#!/usr/bin/env python3
"""
Builds docs/api/code-map.md, the code map, from the pciSeq source.

Two parts. The top is a hand-kept table, MAP below, that goes from a component
of the model (theta, the spatial term, ...) to the functions that compute it. The
rest is generated: every module, class and function under pciSeq/src with its
line and the first line of its docstring, so nothing can be left out. The test
in tests/test_docs_in_sync.py checks that every name in MAP exists and that the
committed page is what this script produces now.

The page is for two readers: a person finding their way round the repo, and the
agent in the viewer chat, which reads the docs first and then the source at the
commit that made the run; with the map it opens the right file at the right line
instead of paging through main.py.

Run it from the website/ folder:  python gen_code_map.py
The deploy workflow runs it too, right before vitepress builds. Standard library
only, it parses the source and does not import the package.
"""

import ast
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent
SRC = REPO / "pciSeq" / "src"
OUT = HERE / "docs" / "api" / "code-map.md"

# Component, a one line description, and the code that computes it as
# (module path under pciSeq/src, "Class.method" or "function"). Order is the
# order of the main loop, then the rest. Keep the names exact, the test checks
# them against the source.
MAP = [
    ("Main loop",
     "One iteration updates the gene counts, rho, eta, theta, gamma, the class probabilities and the spot assignments, in that order. Iterations continue until the stopping rule is met.",
     [("core/main.py", "VarBayes.main_loop"), ("core/main.py", "VarBayes._step"), ("core/main.py", "VarBayes.run")]),
    ("Soft gene counts per cell",
     "Expected counts per cell and gene, each spot contributing its assignment probability. Reported in cellData and used by the class likelihood.",
     [("core/main.py", "VarBayes.geneCount_upd"), ("core/datatypes/cells.py", "Cells.geneCount"),
      ("core/datatypes/spots.py", "Spots.zero_class_counts")]),
    ("rho, the misread density",
     "Per gene rate of spots per unit volume not attributable to any cell. The background candidate in the spot assignment.",
     [("core/main.py", "VarBayes.rho_upd"), ("core/datatypes/genes.py", "Genes.calc_rho"),
      ("core/datatypes/genes.py", "Genes.calc_misread_density"), ("core/datatypes/spots.py", "Spots.misread_density")]),
    ("eta, the gene efficiency",
     "Per gene detection efficiency relative to the Inefficiency setting, scaling the reference expression "
     "to the detected counts. Genes.inefficiency is eta times the Inefficiency setting, not its reciprocal.",
     [("core/main.py", "VarBayes.eta_upd"), ("core/datatypes/genes.py", "Genes.calc_eta"),
      ("core/datatypes/genes.py", "Genes.init_eta"), ("core/datatypes/genes.py", "Genes.inefficiency")]),
    ("theta, the cell scale",
     "Per cell factor scaling the expected counts of a class to the cell's total count, with a Gamma(rTheta, rTheta) prior.",
     [("core/main.py", "VarBayes.theta_upd"), ("core/datatypes/cells.py", "Cells.calc_theta"),
      ("core/datatypes/cells.py", "Cells.init_theta"), ("core/main.py", "VarBayes.init_theta")]),
    ("gamma, the cell-gene scale",
     "Per cell, gene and class factor absorbing overdispersion; integrating it out gives the negative binomial likelihood.",
     [("core/main.py", "VarBayes.gamma_upd"), ("core/datatypes/spots.py", "Spots.init_gamma"),
      ("core/datatypes/spots.py", "Spots.gammaExpectation"), ("core/datatypes/spots.py", "Spots.logGammaExpectation")]),
    ("Expected counts",
     "Expected count of each gene in a cell of a given class: the reference expression times the cell area factor, eta and theta.",
     [("core/utils/likelihood.py", "scaled_exp"), ("core/main.py", "VarBayes.scaled_exp"),
      ("core/datatypes/singleCell.py", "SingleCell.mean_expression")]),
    ("Class probabilities",
     "Per cell, the softmax of the gene log-likelihood, the class prior and the spatial term of each class.",
     [("core/main.py", "VarBayes.cell_to_cellType"), ("core/utils/likelihood.py", "compute_gene_loglikelihood_matrix"),
      ("core/utils/likelihood.py", "negative_binomial_loglikelihood")]),
    ("Class prior",
     "The prior probability of each class, a Dirichlet updated from the current class probabilities.",
     [("core/datatypes/cellClass.py", "CellClass.log_prior"), ("core/datatypes/cellClass.py", "CellClass.ini_prior"),
      ("core/main.py", "VarBayes.dalpha_upd")]),
    ("Spatial term (mrf)",
     "The contribution of the neighbouring cells to a cell's class score, mrf_beta times the weighted class probabilities of the neighbours.",
     [("core/datatypes/cells.py", "Cells.calc_mrf"), ("core/datatypes/cells.py", "Cells.mrf_support"),
      ("core/datatypes/cells.py", "Cells.nearest_neighbours"), ("core/datatypes/cells.py", "Cells.class_pooling")]),
    ("Spot assignments",
     "Per spot, the softmax over the nearby cells and the background of the Gaussian term at the fitted cell centroid and the expression terms.",
     [("core/main.py", "VarBayes.spots_to_cell"), ("core/main.py", "VarBayes.spots_to_cell_numba"),
      ("core/utils/numba_kernels.py", "spots_to_cell_numba_kernel"), ("core/datatypes/spots.py", "Spots.mvn_loglik"),
      ("core/datatypes/spots.py", "Spots.cells_nearby")]),
    ("Stopping rule",
     "Iterations stop when the largest change in a spot assignment probability falls below CellCallTolerance, or at max_iter.",
     [("core/utils/convergence.py", "has_converged"), ("core/utils/convergence.py", "log_iteration_diagnostics")]),
    ("ELBO",
     "The evidence lower bound, computed for diagnostics.",
     [("core/utils/elbo.py", "calc_elbo")]),
    ("Cell geometry",
     "The centre and covariance of each cell's gaussian, a Normal-Inverse-Wishart posterior from the assigned spots. The position term of the spot score is its expected log density.",
     [("core/main.py", "VarBayes.gaussian_upd"), ("core/main.py", "VarBayes.centroid_upd"), ("core/main.py", "VarBayes.cov_upd"),
      ("core/utils/gaussian_model.py", "soft_stats"), ("core/datatypes/spots.py", "Spots.mvn_loglik"), ("core/utils/elbo.py", "geometry_kl")]),
    ("Mean counts per class",
     "Mean count of each gene over the cells assigned to each class, weighted by class probability. Shown in the Gene Expression table of the viewer's cell diagnostics.",
     [("core/datatypes/cells.py", "Cells.mean_gene_reads_per_class"), ("core/datatypes/cells.py", "Cells.gene_reads_per_class")]),
    ("Call diagnostics",
     "Per gene log-likelihood contributions to a cell's class score, and per term scores of a spot's candidate cells.",
     [("core/utils/inspection.py", "check_cell"), ("core/utils/inspection.py", "check_spot"),
      ("core/utils/likelihood.py", "calculate_genes_log_likelihood_contr")]),
    ("Cell label mapping",
     "Conversion between the array rows (1 to nC-1, row 0 the background) and the segmentation labels.",
     [("core/utils/cell_utils.py", "to_internal"), ("core/utils/cell_utils.py", "to_external"),
      ("core/utils/cell_utils.py", "recover_original_labels"), ("core/main.py", "VarBayes.to_internal")]),
    ("Output files",
     "cellData and geneData; the tsv, arrow and spatialdata exports; diagnostics.db.",
     [("core/summary.py", "cells_summary"), ("core/summary.py", "spots_summary"), ("core/io/export.py", "write_data"),
      ("core/io/tsv_export.py", "write_tsv"), ("core/io/arrow_export.py", "write_arrow"),
      ("core/utils/spatialdata_export.py", "write_spatialdata"), ("core/io/diagnostics_db.py", "export_diagnostics")]),
    ("Provenance",
     "pciSeq version, commit and commit date, run date and environment, recorded with every run.",
     [("core/io/provenance.py", "run_metadata"), ("core/io/provenance.py", "serialise")]),
    ("Preprocessing",
     "Conversion of the spots table and the label image into the model's inputs.",
     [("preprocess/main.py", "stage_data")]),
    ("Agent tools",
     "The functions behind the MCP server and the viewer chat. All read diagnostics.db.",
     [("mcp/tools.py", "Run.explain_cell"), ("mcp/tools.py", "Run.explain_spot"), ("mcp/tools.py", "Run.cell_counts"),
      ("mcp/tools.py", "Run.spots_in_cell"), ("mcp/tools.py", "Run.run_info"), ("mcp/tools.py", "narrate_cell"),
      ("mcp/tools.py", "narrate_spot"), ("mcp/tools.py", "open_run")]),
]


def first_line(node):
    doc = ast.get_docstring(node)
    return doc.strip().splitlines()[0].strip() if doc else ""


def modules():
    """Every module under pciSeq/src, sorted, as (relative path, ast tree)."""
    out = []
    for p in sorted(SRC.rglob("*.py")):
        if "__pycache__" in p.parts:
            continue
        out.append((p.relative_to(SRC).as_posix(), ast.parse(p.read_text())))
    return out


def symbols(tree):
    """(qualified name, line, docstring first line) for every class, function and
    method in a module, in source order. Dunders are left out."""
    out = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            out.append((node.name, node.lineno, first_line(node)))
            for sub in node.body:
                if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)) and not sub.name.startswith("__"):
                    out.append(("%s.%s" % (node.name, sub.name), sub.lineno, first_line(sub)))
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("__"):
            out.append((node.name, node.lineno, first_line(node)))
    return out


def index():
    """{module path: {qualified name: line}} for the whole source."""
    return {path: {name: line for name, line, _ in symbols(tree)} for path, tree in modules()}


def missing_map_entries(idx=None):
    """The MAP entries that name something not in the source. Empty when the map
    is current."""
    idx = idx or index()
    return [(path, name) for _, _, refs in MAP for path, name in refs
            if path not in idx or name not in idx[path]]


def cell(text):
    return text.replace("|", "\\|")


def gen_code_map():
    idx = index()
    missing = missing_map_entries(idx)
    if missing:
        raise SystemExit("MAP names not found in the source: %s" % missing)

    out = []
    out.append("---")
    out.append("description: The implementation of each model component in the pciSeq source, and an index of every module and function.")
    out.append("---")
    out.append("")
    out.append("# Code map")
    out.append("")
    out.append("Generated by `website/gen_code_map.py` from the source at the commit this page was built from. "
               "Paths are relative to `pciSeq/src`; line numbers are the first line of each definition.")
    out.append("")
    out.append("## Model components")
    out.append("")
    out.append("| Component | Description | Implementation |")
    out.append("|---|---|---|")
    for what, desc, refs in MAP:
        code = "<br>".join("`%s` `%s` line %d" % (path, name, idx[path][name]) for path, name in refs)
        out.append("| %s | %s | %s |" % (cell(what), cell(desc), code))
    out.append("")
    out.append("## Every module and function")
    out.append("")
    out.append("One block per module, in source order. Methods are listed as `Class.method`; names starting with an underscore are internal.")
    out.append("")
    for path, tree in modules():
        syms = symbols(tree)
        if not syms:
            continue
        head = first_line(tree)
        out.append("<details>")
        out.append("<summary><code>%s</code>%s</summary>" % (path, (" " + cell(head)) if head else ""))
        out.append("")
        out.append("| Name | Line | Description |")
        out.append("|---|---|---|")
        for name, line, doc in syms:
            out.append("| `%s` | %d | %s |" % (name, line, cell(doc)))
        out.append("")
        out.append("</details>")
        out.append("")
    return "\n".join(out)


def main():
    OUT.write_text(gen_code_map())
    print("wrote", OUT.relative_to(REPO))


if __name__ == "__main__":
    main()
