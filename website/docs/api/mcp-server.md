---
description: A program through which an AI assistant reads the files a completed pciSeq run saved.
---

# MCP server

pciSeq includes `pciseq-mcp`, a program through which an AI assistant such as Claude
Code, Gemini CLI or Codex reads the files a completed run saved. The assistant
calls its functions to inspect cell class assignments, spot assignments, counts, scale
factors and images, and answers questions about the run in natural language. MCP
([Model Context Protocol](https://modelcontextprotocol.io)) is the protocol such
assistants use to call outside programs.

It is required only for that means of access; the other three do not use it: an AI
assistant connected to pciSeq Viewer, the chat panel in pciSeq Viewer, and the chat of the
live viewer. They include the same functions. `pciseq-mcp` needs no viewer and no screen,
so it also serves a run on a remote machine. [Accessing the AI agent](../accessing-the-ai-agent.md)
shows the four ways side by side.

The program reads the output folder of the run: `diagnostics.db` and the files for
pciSeq Viewer, which `fit` writes by default, and `cellData.tsv` and `geneData.tsv` when
present. The fitted pickle is not used.

## Quick start

```bash
pip install "pciSeq_3d[mcp] @ git+https://github.com/acycliq/pciSeq_3d.git@dev_3d"
claude mcp add --scope user pciseq -- pciseq-mcp
claude
```

```
Open the run in <output_path> and explain why cell 100 was assigned to its class.
```

These are the three steps described below, in the form for Claude Code.

## Installation

The program is an optional extra. It depends on the `mcp` Python package, which the
rest of pciSeq does not need:

```bash
pip install "pciSeq_3d[mcp] @ git+https://github.com/acycliq/pciSeq_3d.git@dev_3d"
```

This installs the `pciseq-mcp` command, in lowercase. `which pciseq-mcp` prints its
location. Run by hand it prints nothing and waits for an assistant; this is the
expected behaviour, and Ctrl+C exits.

## Registering with an assistant

The program is registered with the assistant once. The assistant then starts it by
itself each time it is needed; it need not be started by hand.

```bash
claude mcp add --scope user pciseq -- pciseq-mcp             # Claude Code
gemini mcp add --scope user pciseq pciseq-mcp                 # Gemini CLI
codex mcp add pciseq -- pciseq-mcp                            # OpenAI Codex CLI and ChatGPT desktop app
```

`pciseq` is the name under which the assistant lists the program and may be any name;
`pciseq-mcp` is the command. `--scope user` makes it available in every folder; without
it Claude Code and Gemini CLI register it for the current folder only. `claude mcp
list`, `gemini mcp list` and `codex mcp list` show what is registered.

### Where the command is found

The assistant runs `pciseq-mcp` as registered, so it must be started from a terminal in
which that command is found. With conda or a virtual environment, that is the
environment pciSeq is installed in. Registering the full path removes the dependence:

```bash
claude mcp add --scope user pciseq -- "$(which pciseq-mcp)"
```

### Assistants without a command

Claude Desktop, Cursor and Zed are configured by adding an entry to a settings file,
with the full path of the command:

```json
{
  "mcpServers": {
    "pciseq": {"command": "/full/path/to/pciseq-mcp", "args": []}
  }
}
```

| Assistant | Settings file | Key |
| --- | --- | --- |
| Claude Desktop | `claude_desktop_config.json`, opened from Settings, Developer, Edit Config | `mcpServers` |
| Cursor | `~/.cursor/mcp.json`, or `.cursor/mcp.json` in the project | `mcpServers` |
| Zed | the settings file | `context_servers` |

Claude Desktop must be restarted afterwards. The commands and locations above are the
assistants' as of October 2026; their own documentation is the reference.

Web applications such as ChatGPT and claude.ai connect to remote servers only and
cannot use a program that runs on the user's machine.

## Usage

Questions are typed to the assistant in natural language. The first call in a session is
`open_run` with the run folder; every other tool refers to the run that is open. With
Claude Code:

```bash
claude
```

```
Open the run in <output_path> and explain why spot 1642419 went to cell 18223.
```

The assistant calls `open_run` on the folder, then `explain_spot`, and answers from the
result. Later questions in the same session need not name the run.

## Tools

*Tool* is the term MCP uses for a function that a program makes available to an
assistant. Each tool has a name, a description and typed arguments, and returns a
result. The assistant reads the descriptions, decides which tools a question calls for,
calls them, and composes its answer from what they return. The tools of `pciseq-mcp`
are listed below.

All tools take and return the cell labels of the input segmentation, the labels shown in
pciSeq Viewer. The internal labels pciSeq assigns when it renumbers a segmentation do
not appear; see [Cell identifiers](./working-with-results.md#cell-identifiers).

### Run

- `open_run(path)`: opens a run. Returns the numbers of cells, spots, genes and classes, the pciSeq version that produced it, and whether containment data is present.
- `run_info()`: provenance and settings. The pciSeq version, commit and commit date, the run date, the Python and package versions, the mean cell radius, the resolved settings, the number of iterations, and whether the loop converged.
- `metadata(key=None)`: the keys of the `metadata` table of `diagnostics.db`, or the value of one, for quantities without a tool of their own such as the cell type definitions.

### Cells

- `cell(label)`: the class probabilities, top genes, total counts and scale factor of one cell, and the neighbours that enter its spatial term.
- `cell_shape(label)`: the fitted shape of the cell's Gaussian: round or elongated, the radii along its own axes, the direction of the longest, each radius against the mean cell radius, and the covariance. Runs saved before October 2026 do not keep the covariance.
- `explain_cell(label, vs_class=None, top_n=10)`: the score of the assigned class against a second class, split into the gene log-likelihood, the class prior and the spatial term, with the genes contributing most to each side. Each gene is reported with the cell's count and the mean count over the cells assigned to either class. The shared genes are the cell's largest counts that the two classes fit about equally. `vs_class` defaults to the runner-up.
- `cell_counts(label, gene=None)`: the counts of a cell, in total or for one gene.
- `spots_in_cell(label, gene=None)`: the spots whose pixel lies inside the cell's segmentation mask.
- `spots_of_cell(label, min_prob=None, gene=None)`: the spots whose most probable parent is the cell, each with its probability. With `min_prob`, every spot with probability above it; with `gene`, only that gene's spots.
- `cell_row(label)`: the `cellData.tsv` row of one cell.
- `theta(label)`: the cell scale factor, overall and under each of the top classes.
- `gamma(label, gene=None)`: the cell-gene scale factors under the assigned class, for every gene or for one.
- `neighbours(label)`: the cells entering the spatial term of the cell, with their classes and, when the files for pciSeq Viewer are present, the distances between their fitted centroids.
- `class_counts(min_counts=None)`: the number of cells per class, hard and soft.
- `class_prior()`: the prior weight of every class at the end of the run, the probability that a cell is of the class before its genes are considered. Zero is listed first, the rest by weight.
- `find_cells(class_name=None, plane=None, min_counts=None, top_two_within=None, n=50)`: the cells matching the filters, by assigned class, plane, total counts, or the gap between the top two classes.

### Spots

- `spot(spot_id)`: the gene, position, plane and candidate cells of one spot, with their probabilities.
- `explain_spot(spot_id)`: one row per candidate cell plus the background: the six score terms, their sum and the resulting probability.
- `spot_row(spot_id)`: the `geneData.tsv` row of one spot.

### Genes

- `gene(name)`: for one gene, eta (its detection efficiency relative to the Inefficiency setting) and eta times the Inefficiency setting; its misread density, its spot totals, its soft counts per class, and the cells holding most of it.

### Images

- `cell_image(label, context=False, plane=None, width=1200, channel=None, neighbours=False, save_as=None, mbtiles=None)`: a rendering of the cell on the tissue image, stitched from the `.mbtiles` of pciSeq Viewer: a close-up with the cell outlined in red and the other cells of the plane in blue, or with `context=True` the whole plane with the cell marked. `neighbours=True` outlines only the cells that enter the spatial term. `save_as` writes the PNG to a file.
- `plane_image(plane=None, bbox=None, width=1200, channel=None, save_as=None, mbtiles=None)`: the tissue image of one plane, whole or restricted to `bbox` in image pixels. The plane defaults to the middle of the stack. When the run has several background images, `channel` selects one by name; without it the tool lists them.

### Arithmetic

- `calculate(expression)`: evaluates an arithmetic expression: numbers, `+ - * / **`, brackets, and `exp`, `log`, `log10`, `sqrt`, `abs` and `round`. Anything else is refused. The assistant is instructed to use it rather than compute in its reply, and not to use it to form quantities the other tools do not define.

### Source code

- `list_source(dir='')`: a folder of the pciSeq source, folders and python files, at the commit that made the open run.
- `read_source(path, start_line=1)`: a file of the source with line numbers, at most 400 lines from `start_line`, at the commit that made the open run. A path as the repository has it (`pciSeq/src/core/main.py`) or as the [code map](./code-map.md) writes it (`core/main.py`).

Both serve the code that produced the run's numbers, read from GitHub at the run's commit, which needs a network connection. A run from before the commit was recorded is refused.

### Documentation

- `docs(query='', page='', n=5)`: the paragraphs of this documentation matching a keyword query, each with its page and heading; with `page`, one whole page; with neither, the list of pages.

Every page is also available as a resource named `pciseq-docs://<page>`, and
`pciseq-docs://index` lists them.

With a run open, the pages are those saved inside its `diagnostics.db` when it was
fitted, so that a run is described by the documentation of the version that produced
it. A run that carries no pages is refused: `docs` and the resources report that no
documentation is available. With no run open, the pages are those of the installed
pciSeq, packed into the wheel at build time from `website/docs`; in a repository
checkout the folder itself is read.

`explain_cell` and `explain_spot` return the same numbers as
[`check_cell`](./reference.md#check-cell) and [`check_spot`](./reference.md#check-spot),
which [Explaining the calls](../explaining-the-calls/overview.md) works through by hand.

## Return values

All tools return JSON-serialisable dictionaries, except `cell_image` and `plane_image`,
which return the PNG image and a JSON description of it. Errors, such as an unknown label or a
run without containment data, are returned as tool errors with a message rather than
raised.

`explain_cell` and `explain_spot` include a `narrative` field: the result in prose,
generated from the numbers, with log-likelihood differences expressed as odds. The
narrative is deterministic and is the same in the chat panel of pciSeq Viewer.

Two tools count the spots of a cell and give different numbers by design. `cell_counts`
returns the sum of the assignment probabilities of the spots, the quantity stored as
`CellGeneCount` in `cellData`; it is not an integer. `spots_in_cell` returns the number
of spots whose pixel lies inside the cell's mask; no probability enters it. The model
scores a spot's position under the cell's Gaussian, centred on the fitted cell centroid, not against the mask, so a spot outside
every mask can still be assigned to a cell. Each result states which quantity it is.

Probabilities in `cellData.tsv` and `geneData.tsv` are stored to three decimals and
`cell_row` and `spot_row` return them as stored. `explain_spot` recomputes them at full
precision.

## Examples

The run is the Espio section used throughout the documentation.

**Why is cell 2413 `048 RHP-COA Ndnf Gaba`?** The assistant calls `explain_cell(2413)`.
The result:

```
                    assigned   compared
gene log-likelihood  -674.94    -708.54
class prior            -4.33      -4.33
spatial term            4.80       0.01

for the assigned class   Ndnf +15.3   Rgs5 +15.3   Ttr +10.7   Vip +4.9
for the compared class   Npy   -7.6   Kit  -7.2   Rgs10 -3.9
```

and its `narrative`:

> Cell 2413 was called 048 RHP-COA Ndnf Gaba, with probability 1.00. The closest
> alternative was 047 Sncg Gaba, at less than 0.01. pciSeq decides a cell's class from
> three things: how well its gene counts match what each class expresses according to
> the cell type definitions (the gene log-likelihood), how common each class is to begin
> with (the prior), and what the neighbouring cells were called (the spatial term). The
> class with the best total receives the highest probability. The genes point to 048
> RHP-COA Ndnf Gaba,
> overwhelmingly, beyond any doubt. The strongest evidence comes from Ndnf, Rgs5 and
> Ttr: the cell holds these in amounts that fit what the cell type definitions give for
> a 048 RHP-COA Ndnf Gaba cell, and not for a 047 Sncg Gaba cell. A few genes, Npy, Kit
> and Rgs10, look more like 047 Sncg Gaba, but they are outweighed. The prior treats the
> two classes alike. The neighbouring cells are mostly 048 RHP-COA Ndnf Gaba, which
> strengthens the call. So the genes settled it, and the neighbourhood agreed.

When the spatial term outweighs the genes, the narrative states it. Cell 4308 is
assigned `022 L5 ET CTX Glut` although its genes favour `006 L4/5 IT CTX Glut` by about
60 to one; the neighbouring cells are L5 ET, and the spatial term determines the assignment.

**Why did spot 1642419 go to cell 18223?** The assistant calls `explain_spot(1642419)`.
The top three candidates:

```
cell         class          Gaussian fit    class expr.  cell scale   cell-gene   prob
18223        037 DG Glut          -10.41           0.94       -0.66        0.18   0.74
21574        037 DG Glut          -12.81           0.94        0.10        0.07   0.13
17371        037 DG Glut          -12.07           0.94       -0.30       -0.54   0.10
background                                                                        0.01   misread -14.92
```

and its `narrative`:

> Spot 1642419 is a Synpr spot. It was assigned to cell 18223, fairly confidently, with
> probability 0.74. The next candidates are cell 21574 (0.13), cell 17371 (0.10), and
> the chance it is a misread is 0.01. pciSeq weighs each nearby cell on two things: how
> close the spot is to the cell's centroid (the Gaussian fit), and how well a Synpr spot
> fits that cell, which combines whether the cell's class expresses Synpr (class
> expression), whether the cell holds more reads overall than its class predicts (cell
> scale), and whether it already holds more Synpr than its class predicts (cell-gene
> scale). The background is scored on how often Synpr spots turn out to be misreads. The
> candidate with the best total receives the highest probability. Cell 18223 is the
> nearest candidate: about 5 to one over cell 17371
> and about 11 to one over cell 21574 on position alone. Every candidate is a 037 DG
> Glut cell, and on class alone Synpr fits them about equally; what separates them is
> how much each already holds. Between the top two, cell 21574 holds more reads overall
> than its class predicts. So cell 21574 is the better fit for the gene, slightly, but cell
> 18223 is closer, about 11 to one, and distance decides it.

The probabilities match [Table 3.1](../explaining-the-calls/why-a-spot-got-its-cell.md#table-3-1).
All candidates are `037 DG Glut` and their class expression terms are equal (0.94), so
that term does not separate them; cell 18223
is the nearest and the Gaussian term determines the assignment. The narrative states this as *distance
decides it*; when a farther cell is assigned the spot because its class expresses the gene, it
states *expression decides it*.

**How many Ndnf counts does cell 2413 have?** `cell_counts(2413, gene='Ndnf')` returns
7.77, marked as a soft count.

**Which spots belong to cell 18223?** `spots_of_cell(18223)` returns the 31 spots whose
most probable parent is the cell, with probabilities from 0.810 to 0.295. `cell_row(18223)`
lists 506 spots under `spot_id`, every spot with probability above 0.0001 on the cell;
their probabilities sum, to within that cut-off, to the cell's 38.5 counts.
