#!/usr/bin/env python3
"""
Builds the API reference markdown for the VitePress docs straight from the
pciSeq source, so the docs don't drift away from the code.

Two things get generated into docs/api/:
  - reference.md      the public functions (fit, cell_type, stage_data, ...),
                      pulled from their signatures and docstrings.
  - configuration.md  the opts/config dictionary, pulled from the comments
                      sitting above each key in pciSeq/config.py.

Run it from the website/ folder:  python gen_api.py
The deploy workflow runs it too, right before vitepress builds.
"""

import ast
import inspect
import json
import pathlib
import re
import textwrap

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent
OUT = HERE / "docs" / "api"

# Extra things to document that discover_public_api() does not catch on its own.
# Each entry is (module path, attribute name).
#   - stage_image is reachable as pciSeq.stage_image, but only through a
#     conditional import (it depends on libvips), so the top-level scan misses it.
#   - VarBayes is the model object that cell_type/fit build and return.
# stage_image is a plain function so it goes with the others; VarBayes (a class)
# stays last.
EXTRA_API = [
    ("pciSeq.src.tiling.stage_image", "stage_image"),
    ("pciSeq.src.core.main", "VarBayes"),
]

# numpydoc sections that list "name : type" entries we want as bullet lists.
PARAM_SECTIONS = {"Parameters", "Returns", "Raises", "Yields", "Attributes"}

# Classes have a lot of internal methods we don't want in the docs. For these
# classes, list exactly which members to show (in order). Entries are either a
# plain method name on the class itself, or a member reached through one of the
# instance's attributes, written as "attr.method" with where it actually lives.
CLASS_MEMBERS = {
    "VarBayes": {
        # a short note rendered right after the class summary, so readers know
        # how to actually get hold of a fitted model before the attributes and
        # methods below make sense.
        "note": (
            "::: tip Obtaining a fitted instance\n"
            "`VarBayes` is not instantiated directly in normal use. "
            "[`fit`](#fit) and [`cell_type`](#cell-type) construct and run it. "
            "`cell_type` returns the fitted instance; `fit` does not, but when "
            "`save_data=True` (the default) the fitted model is serialised to "
            "`<output_path>/data/debug/pciSeq.pickle` (`output_path` defaults to a "
            "temporary directory). The attributes and methods below operate on a "
            "loaded instance; [Working with results](./working-with-results) walks "
            "through the main ones with examples.\n\n"
            "```python\n"
            "import pandas as pd\n\n"
            "obj = pd.read_pickle('<output_path>/data/debug/pciSeq.pickle')\n\n"
            "obj.metadata\n"
            "obj.check_cell(my_label=42, user_class='Astro')\n"
            "```\n"
            ":::"
        ),
        # plain instance attributes worth surfacing. Attributes have no
        # docstrings, so the description is written here by hand.
        "attributes": [
            {
                "name": "metadata",
                "type": "dict",
                "desc": "Provenance recorded when the model is built, saved "
                        "alongside the results so a run can be traced back to the "
                        "code that produced it. Contains:",
                "fields": [
                    ("version", "the pciSeq version"),
                    ("branch", "the git branch"),
                    ("commit", "the git commit hash"),
                    ("build_date", "the package build date"),
                    ("created_at", "a UTC timestamp for when the run was created"),
                ],
            },
        ],
        "methods": [
            "check_spot",
            "check_cell",
            "read_tsv",
            "heatmap_counts_per_class",
        ],
        # reached as varBayes.cells.<method>, but defined on the Cells class.
        # gene_reads_per_class is the numerator, so it goes first.
        "attr_methods": [
            {
                "display": "cells.gene_reads_per_class",
                "module": "pciSeq.src.core.datatypes.cells",
                "owner": "Cells",
                "method": "gene_reads_per_class",
            },
            {
                "display": "cells.mean_gene_reads_per_class",
                "module": "pciSeq.src.core.datatypes.cells",
                "owner": "Cells",
                "method": "mean_gene_reads_per_class",
            },
        ],
    },
}


def discover_public_api():
    """Read pciSeq/__init__.py and return the names it pulls into the package
    namespace, i.e. everything you can reach as pciSeq.<name>.

    We parse the file instead of importing it, so adding a new export there
    automatically shows up in the docs with no edit here.
    """
    init = (REPO / "pciSeq" / "__init__.py").read_text()
    found = []
    for node in ast.parse(init).body:
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("pciSeq"):
            for alias in node.names:
                name = alias.asname or alias.name
                if name.startswith("_"):  # skip __version__ and friends
                    continue
                found.append((node.module, name))
    return found


def _dedent_body(lines):
    """Drop the leading blank lines and dedent a block of docstring lines."""
    text = "\n".join(lines)
    return textwrap.dedent(text).strip("\n")


def parse_numpydoc(doc):
    """Split a NumPy-style docstring into (summary, [(header, body_lines)]).

    A section starts with a title line followed by a line of dashes, e.g.

        Parameters
        ----------
    """
    lines = (doc or "").expandtabs().splitlines()
    summary, sections = [], []
    i, n = 0, len(lines)
    current_header, current_body = None, []

    def flush():
        if current_header is not None:
            sections.append((current_header, current_body[:]))

    while i < n:
        line = lines[i]
        nxt = lines[i + 1] if i + 1 < n else ""
        is_header = (
            line.strip()
            and set(nxt.strip()) == {"-"}
            and len(nxt.strip()) >= len(line.strip())
        )
        if is_header:
            flush()
            current_header = line.strip()
            current_body = []
            i += 2
            continue
        if current_header is None:
            summary.append(line)
        else:
            current_body.append(line)
        i += 1
    flush()
    return _dedent_body(summary), sections


def render_param_section(body_lines):
    """Render a Parameters/Returns/... block as a markdown bullet list.

    Entries look like:
        name : type
            description, possibly several lines
    """
    body = _dedent_body(body_lines)
    out = []
    entry_desc = []

    def flush_desc():
        if entry_desc:
            desc = " ".join(d.strip() for d in entry_desc if d.strip())
            if desc:
                out.append(f"  {desc}")
            entry_desc.clear()

    for line in body.splitlines():
        # a new entry starts at column 0 and is not indented continuation text
        if line and not line[0].isspace():
            flush_desc()
            if " : " in line:
                name, _, typ = line.partition(" : ")
                out.append(f"- **`{name.strip()}`** *({typ.strip()})*")
            else:
                out.append(f"- **`{line.strip()}`**")
        else:
            entry_desc.append(line)
    flush_desc()
    return "\n".join(out)


_MODULE_CACHE = {}


def _module(module_path):
    """Read and parse a module's source once, return (source_text, tree).

    We read the source straight off disk and parse it, instead of importing
    pciSeq. That keeps this script dependency-free so CI can run it without
    installing numpy/pandas/scipy.
    """
    if module_path not in _MODULE_CACHE:
        rel = pathlib.Path(*module_path.split(".")).with_suffix(".py")
        src = (REPO / rel).read_text()
        _MODULE_CACHE[module_path] = (src, ast.parse(src))
    return _MODULE_CACHE[module_path]


def _find_def(module_path, name):
    """Find a top-level function or class, return (node, module_source)."""
    src, tree = _module(module_path)
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.name == name:
            return node, src
    raise RuntimeError(f"could not find {name} in {module_path}")


def _raw_docstring(node, source):
    """The docstring as written in the source, backslashes intact.

    ast.get_docstring evaluates the string, so a docstring that isn't a raw
    string turns \\frac into a form-feed and loses the backslash. We pull the
    literal source text instead, then dedent it like get_docstring would.
    """
    body = getattr(node, "body", [])
    if not (body and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)):
        return None
    seg = ast.get_source_segment(source, body[0].value)
    if seg is None:
        return ast.get_docstring(node)
    seg = seg.strip()
    seg = re.sub(r'^[rRbBuUfF]*("""|\'\'\'|"|\')', "", seg)
    seg = re.sub(r'("""|\'\'\'|"|\')$', "", seg)
    # cleandoc handles the first line having no indent, same as get_docstring
    return inspect.cleandoc(seg)


def _signature(node, display_name=None):
    """Build a readable 'name(args) -> return' string from a function node.

    ast.unparse keeps the annotations exactly as written in the source
    (pd.DataFrame, Tuple[...]), which reads better than the fully-qualified
    names inspect would give us.
    """
    args = ast.unparse(node.args)
    sig = f"{display_name or node.name}({args})"
    if node.returns is not None:
        sig += f" -> {ast.unparse(node.returns)}"
    return sig


def _is_property(node):
    return any(
        (isinstance(d, ast.Name) and d.id == "property")
        for d in node.decorator_list
    )


def _find_method(class_node, method_name):
    for n in class_node.body:
        if isinstance(n, ast.FunctionDef) and n.name == method_name:
            return n
    return None


def _fix_math_delims(text):
    """Some docstrings use LaTeX \\( \\) and \\[ \\] delimiters. VitePress math
    wants $ and $$. We also strip the spaces just inside the delimiters, because
    inline math like '$ w $' (space right after the $) does not render."""
    text = re.sub(r"\\\(\s*(.+?)\s*\\\)", r"$\1$", text, flags=re.DOTALL)
    text = re.sub(r"\\\[\s*(.+?)\s*\\\]", r"$$\1$$", text, flags=re.DOTALL)
    return text


def _render_method(node, display_name, source):
    """A single #### method block: signature plus its docstring."""
    sig = _signature(node, display_name=display_name)
    sig = sig.replace("(self, ", "(").replace("(self)", "()")
    tag = " *(property)*" if _is_property(node) else ""
    parts = [f"#### `{display_name}`{tag}", "", "```python", sig, "```", ""]
    msum, msec = parse_numpydoc(_raw_docstring(node, source))
    parts += _render_docbody(msum, msec)
    return parts


def _render_docbody(summary, sections):
    """The shared body: summary plus each docstring section. Section labels are
    bold rather than headings so they don't clutter the page outline."""
    parts = []
    if summary:
        parts += [_fix_math_delims(summary), ""]
    for header, body in sections:
        parts.append(f"**{header}**")
        parts.append("")
        if header in PARAM_SECTIONS:
            parts.append(render_param_section(body))
        else:
            parts.append(_fix_math_delims(_dedent_body(body)))
        parts.append("")
    return parts


def render_function(node, heading, qualpath, source):
    parts = [heading, "", f"`{qualpath}`", "", "```python", _signature(node), "```", ""]
    summary, sections = parse_numpydoc(_raw_docstring(node, source))
    parts += _render_docbody(summary, sections)
    return "\n".join(parts).rstrip() + "\n"


def render_class(node, name, qualpath, source):
    # constructor signature, shown under the class name itself
    init = next(
        (n for n in node.body if isinstance(n, ast.FunctionDef) and n.name == "__init__"),
        None,
    )
    ctor = _signature(init, display_name=name) if init else f"{name}(...)"
    # drop the leading 'self, ' and the '-> None' a constructor never returns
    ctor = ctor.replace("(self, ", "(").replace("(self)", "()")
    ctor = ctor.split(" -> ")[0]

    parts = [f"## `{name}`", "", f"`{qualpath}`", "", "```python", ctor, "```", ""]
    summary, sections = parse_numpydoc(_raw_docstring(node, source))
    parts += _render_docbody(summary, sections)

    spec = CLASS_MEMBERS.get(name)

    # a short orientation note (e.g. how to get hold of a fitted instance)
    if spec and spec.get("note"):
        parts += [spec["note"], ""]

    # hand-written attributes (plain instance attributes have no docstrings)
    if spec and spec.get("attributes"):
        parts += ["### Attributes", ""]
        for attr in spec["attributes"]:
            # bullets get no auto anchor, so add one for the sidebar link to hit
            anchor = f'<a id="{_slug(attr["name"])}"></a>'
            parts.append(f"- {anchor}**`{attr['name']}`** *({attr['type']})*")
            parts.append(f"  {attr['desc']}")
            # the dict keys (or similar) as a nested list, so they aren't crammed
            for key, meaning in attr.get("fields", []):
                parts.append(f"    - `{key}`: {meaning}")
        parts.append("")

    # which methods to show: a hand-picked list for noisy classes, otherwise
    # every public method.
    method_blocks = []
    if spec:
        for mname in spec.get("methods", []):
            m = _find_method(node, mname)
            if m is None:
                raise RuntimeError(f"{name}.{mname} listed in CLASS_MEMBERS but not found")
            method_blocks += _render_method(m, mname, source)
        for extra in spec.get("attr_methods", []):
            owner, owner_src = _find_def(extra["module"], extra["owner"])
            m = _find_method(owner, extra["method"])
            if m is None:
                raise RuntimeError(f"{extra['owner']}.{extra['method']} not found")
            method_blocks += _render_method(m, extra["display"], owner_src)
    else:
        for m in node.body:
            if isinstance(m, ast.FunctionDef) and not m.name.startswith("_"):
                method_blocks += _render_method(m, m.name, source)

    if method_blocks:
        parts += ["### Methods", ""] + method_blocks
    return "\n".join(parts).rstrip() + "\n"


def render_member(module_path, name):
    """Render whichever it is, a function or a class."""
    node, source = _find_def(module_path, name)
    qualpath = f"{module_path}.{name}"
    if isinstance(node, ast.ClassDef):
        return render_class(node, name, qualpath, source)
    return render_function(node, f"## `{name}`", qualpath, source)


def _slug(name):
    """Match VitePress's heading-id slugify so sidebar anchors line up: lowercase,
    and turn every run of punctuation/space (dots, underscores, etc.) into a
    single hyphen. This is what gives `cells.mean_gene_reads_per_class` the id
    `cells-mean-gene-reads-per-class`, same as VitePress."""
    s = name.strip().lower()
    s = re.sub(r"[\s~`!@#$%^&*()\-_+=\[\]{}|\\;:\"'<>,.?/]+", "-", s)
    return re.sub(r"-{2,}", "-", s).strip("-")


def _class_subitems(spec):
    """Sidebar sub-items for a curated class: its attributes, then methods, then
    attribute-reached methods, each pointing at its anchor on the reference page."""
    items = []
    for attr in spec.get("attributes", []):
        items.append({"text": attr["name"], "link": f"/api/reference#{_slug(attr['name'])}"})
    for mname in spec.get("methods", []):
        items.append({"text": mname, "link": f"/api/reference#{_slug(mname)}"})
    for extra in spec.get("attr_methods", []):
        items.append({"text": extra["display"], "link": f"/api/reference#{_slug(extra['display'])}"})
    return items


def gen_nav():
    """The list of API members for the sidebar, so the 'Functions' group can
    expand into one entry per member. Curated classes (e.g. VarBayes) become a
    nested group of their attributes and methods. config.mts imports this json."""
    members = discover_public_api() + EXTRA_API
    nav = []
    for _module_path, name in members:
        entry = {"text": name, "link": f"/api/reference#{_slug(name)}"}
        spec = CLASS_MEMBERS.get(name)
        if spec:
            entry["collapsed"] = True
            entry["items"] = _class_subitems(spec)
        nav.append(entry)
    return nav


def gen_reference():
    members = discover_public_api() + EXTRA_API
    chunks = [
        "# API reference",
        "",
        "::: warning Auto-generated",
        "This page is generated from the pciSeq source by `website/gen_api.py`.",
        "Edit the docstrings in the source, not this file.",
        ":::",
        "",
        "Everything here is reachable as `pciSeq.<name>` (plus `VarBayes`, the",
        "model object that [`fit`](#fit) and [`cell_type`](#cell-type) build and",
        "return). The main entry point is [`fit`](#fit). For what the output",
        "DataFrames hold and worked examples of the model attributes, see",
        "[Working with results](./working-with-results).",
        "",
    ]
    for module_path, name in members:
        chunks.append(render_member(module_path, name))
        chunks.append("")
    return "\n".join(chunks).rstrip() + "\n"


def _comment_block_above(src_lines, key_lineno):
    """Grab the run of `# ...` comment lines sitting right above a config key.

    key_lineno is 1-based (ast convention). Walk upward over comment lines,
    stop at the first blank line or code line.
    """
    comments = []
    idx = key_lineno - 2  # line directly above the key, 0-based
    while idx >= 0:
        stripped = src_lines[idx].strip()
        if stripped.startswith("#"):
            comments.append(stripped.lstrip("#").strip())
            idx -= 1
        else:
            break
    comments.reverse()
    return comments


def gen_configuration():
    cfg_path = REPO / "pciSeq" / "config.py"
    src = cfg_path.read_text()
    src_lines = src.splitlines()
    tree = ast.parse(src)

    # find the DEFAULT = {...} assignment
    default_dict = None
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and tgt.id == "DEFAULT":
                    default_dict = node.value
    if not isinstance(default_dict, ast.Dict):
        raise RuntimeError("could not find the DEFAULT dict in config.py")

    chunks = [
        "# Configuration (opts)",
        "",
        "::: warning Auto-generated",
        "This page is generated from the comments in `pciSeq/config.py` by",
        "`website/gen_api.py`. Edit the comments in that file, not this page.",
        ":::",
        "",
        "Pass any of these as an `opts` dictionary to [`fit`](./reference#fit).",
        "Anything you leave out keeps its default shown below.",
        "",
        "```python",
        "import pciSeq",
        "opts = {'max_iter': 500, 'CellCallTolerance': 0.01}",
        "cellData, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=ref, opts=opts)",
        "```",
        "",
    ]

    for key_node, val_node in zip(default_dict.keys, default_dict.values):
        if not isinstance(key_node, ast.Constant):
            continue
        key = key_node.value
        default_repr = ast.unparse(val_node)
        comments = _comment_block_above(src_lines, key_node.lineno)
        chunks.append(f"### `{key}`")
        chunks.append("")
        chunks.append(f"**Default:** `{default_repr}`")
        chunks.append("")
        if comments:
            chunks.append("\n".join(comments))
            chunks.append("")
    return "\n".join(chunks).rstrip() + "\n"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "reference.md").write_text(gen_reference())
    (OUT / "configuration.md").write_text(gen_configuration())
    nav_path = HERE / "docs" / ".vitepress" / "api-nav.json"
    nav_path.write_text(json.dumps(gen_nav(), indent=2) + "\n")
    print(f"wrote {OUT/'reference.md'}")
    print(f"wrote {OUT/'configuration.md'}")
    print(f"wrote {nav_path}")


if __name__ == "__main__":
    main()
