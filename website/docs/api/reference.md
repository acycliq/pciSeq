# API reference

::: warning Auto-generated
This page is generated from the pciSeq source by `website/gen_api.py`.
Edit the docstrings in the source, not this file.
:::

Everything here is reachable as `pciSeq.<name>` (plus `VarBayes`, the
model object that [`fit`](#fit) and [`cell_type`](#cell-type) build and
return). The main entry point is [`fit`](#fit).

## `fit`

`pciSeq.app.fit`

```python
fit(*args, **kwargs) -> Tuple[pd.DataFrame, pd.DataFrame]
```

Main entry point for pciSeq cell typing analysis.

You can pass spots and coo either as the first two positional arguments or as
keywords. The keywords below are the preferred way since they read better.

**Parameters**

- **`spots`** *(pd.DataFrame)*
  The spots to assign. Needs the columns 'gene_name', 'x' and 'y', plus 'z_plane' for 3D data.
- **`coo`** *(list of scipy.sparse.coo_matrix)*
  The label image, one sparse matrix per z-plane. A list with more than one plane is treated as 3D.
- **`scRNAseq`** *(pd.DataFrame, optional)*
  Single-cell reference data used to annotate the cell types. Leave it out to run without a reference.
- **`opts`** *(dict, optional)*
  Any config values you want to override, e.g. {'max_iter': 500}. See the configuration page for the full list of keys and their defaults.

**Returns**

- **`cellData`** *(pd.DataFrame)*
  Cell typing results and metadata, one row per cell.
- **`geneData`** *(pd.DataFrame)*
  Gene assignment results, one row per spot.

**Raises**

- **`ValueError`**
  If spots or coo are missing or invalid.
- **`RuntimeError`**
  If the cell typing algorithm fails. Non-convergence on its own only logs a warning, it does not raise.


## `cell_type`

`pciSeq.app.cell_type`

```python
cell_type(cells: pd.DataFrame, spots: pd.DataFrame, scRNAseq: Optional[pd.DataFrame], config: Dict[str, Any], viewer: Optional[Any]=None) -> Tuple[pd.DataFrame, pd.DataFrame, VarBayes]
```

Perform cell typing using Variational Bayes algorithm.

**Parameters**

- **`cells`** *(pd.DataFrame)*
  Preprocessed cell data containing cell locations and boundaries
- **`spots`** *(pd.DataFrame)*
  Preprocessed spot data containing gene expressions and coordinates
- **`scRNAseq`** *(Optional[pd.DataFrame])*
  Single-cell RNA sequencing reference data. Can be None if not using reference data
- **`config`** *(Dict[str, Any])*
  Configuration dictionary containing algorithm parameters
- **`viewer`** *(optional)*
  A running RealtimeViewerServer to stream iterations to, or None. When given, it is bound to the model via viewer.attach(varBayes).

**Returns**

- **`Tuple[pd.DataFrame, pd.DataFrame, VarBayes]`**
  - cellData: DataFrame containing cell typing results - geneData: DataFrame containing gene assignment results - varBayes: The fitted VarBayes model instance

**Raises**

- **`ValueError`**
  If input data is invalid or incompatible
- **`RuntimeError`**
  If the cell typing algorithm fails (non-convergence only logs a warning)


## `stage_data`

`pciSeq.src.preprocess.main.stage_data`

```python
stage_data(spots: pd.DataFrame, coo: List[coo_matrix], cfg: Dict) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
```

Process spots and label images for cell typing analysis.

**Parameters**

- **`spots`** *(pd.DataFrame)*
  Spot data with columns: ['gene_name', 'x', 'y', 'z_plane']
- **`coo`** *(List[coo_matrix])*
  List of sparse matrices containing cell segmentation
- **`cfg`** *(Dict)*
  Configuration dictionary with processing parameters

**Returns**

- **`cells`** *(pd.DataFrame)*
  Cell properties including position and size
- **`borders_future`** *(Future)*
  Future resolving to (cell_boundaries, cell_boundaries_list). Border extraction runs in the background and only blocks when .result() is called.
- **`processed_spots`** *(pd.DataFrame)*
  Processed spots with cell assignments

**Note**

The label remapping (label_map) and image dimensions (img_dim) are written
into cfg as runtime state, the same way Config.set_runtime_attrs adds is3D.


## `attach_to_log`

`pciSeq.src.core.logger.attach_to_log`

```python
attach_to_log()
```

exists only for backwards compatibility.
Replaced by setup_logger


## `setup_logger`

`pciSeq.src.core.logger.setup_logger`

```python
setup_logger(level=None)
```

Configure pciSeq logging with colored console output.

WARNING: This function clears all existing root logger handlers and replaces
them with pciSeq's own handler. If pciSeq is embedded inside a larger
application that has its own logging setup, do NOT call this function.
The parent application's handlers will be wiped out. Only call setup_logger()
when pciSeq is the top-level application.

Args:
    level: logging level (e.g. logging.DEBUG, logging.INFO). Defaults to INFO.

Returns:
    The configured 'pciSeq' logger instance.


## `VarBayes`

`pciSeq.src.core.main.VarBayes`

```python
VarBayes(cells_df: pd.DataFrame, spots_df: pd.DataFrame, scRNAseq: pd.DataFrame, config: Dict[str, Any])
```

Implements Variational Bayes algorithm for spatial transcriptomics analysis.

This class performs cell type assignment and spot-to-cell mapping using a
probabilistic model with variational inference.

Args:
    cells_df: DataFrame containing cell information
    spots_df: DataFrame containing spot information
    scRNAseq: Single-cell RNA sequencing reference data
    config: Configuration dictionary containing algorithm parameters

### Methods

#### `check_spot`

```python
check_spot(spot_id)
```

Break down the spot-to-cell score for a single spot.

For one spot, shows how its assignment score splits across the candidate
cells (plus the background/misread option). The score is the sum of the
location term, the attention, the expression fluctuation, the cell and gene
efficiency terms and the inside-cell bonus. It draws the score and
probability charts and hands back the same numbers as a table.

**Parameters**

- **`spot_id`** *(int)*
  The spot id (its index in the spots table) to look at.

**Returns**

- **`pd.DataFrame`**
  One row per candidate cell plus a background row, with each score term, the misread value, and their sum.

#### `check_cell`

```python
check_cell(my_label, user_class, top_n=10, show_plot=True)
```

Compare two cell types for one cell, gene by gene.

Takes a cell and digs into why the model typed it the way it did, by
putting the class it was assigned head to head against a class you pick.
It pulls the per-gene log-likelihood contributions for both classes so you
can see which genes pushed the call one way or the other, and can draw the
comparison for you.

**Parameters**

- **`my_label`** *(int)*
  The cell label (cell number) to look at. If the segmentation labels were renumbered internally you still pass the original label here, it gets mapped for you.
- **`user_class`** *(str)*
  The other class you want to weigh the assigned class against.
- **`top_n`** *(int, optional)*
  How many genes to show at each end, i.e. the genes that argue hardest for the assigned class and the ones that argue hardest for your class. Default is 10.
- **`show_plot`** *(bool, optional)*
  If True (the default) draw the comparison figure. Set it to False if you just want the tables back.

**Returns**

- **`gene_expression_data`** *(pd.DataFrame)*
  One row per selected gene. Columns hold the mean counts across all cells of each class, the negative-binomial counts the model expects for this cell under each class, and this cell's own observed counts.
- **`contributions`** *(pd.DataFrame)*
  Per-gene log-likelihood for the two classes, so you can read off how much each gene pulled.
- **`fig`** *(matplotlib.figure.Figure or None)*
  The comparison figure, or None when show_plot is False.

#### `read_tsv`

```python
read_tsv(filepath)
```

Read a tsv file that pciSeq wrote back into a DataFrame.

Same as a tab-separated pandas.read_csv, but it also turns the columns
that pciSeq saved as text (lists, dicts or tuples written out as strings
like "[1, 2, 3]") back into real Python objects.

**Parameters**

- **`filepath`** *(str)*
  Path to the tsv file, e.g. a cellData.tsv or geneData.tsv that a run produced.

**Returns**

- **`pd.DataFrame`**
  The file contents, with the list/dict/tuple columns parsed back from their string form.

#### `heatmap_counts_per_class`

```python
heatmap_counts_per_class()
```

Draw the mean-gene-reads-per-class heatmap.

Builds an interactive Plotly heatmap of the average reads for every gene
in every cell class (genes down the rows, classes across the columns) and
shows it. Handy for eyeballing which genes mark which classes. This one is
display only, it does not return anything.

**Returns**

- **`None`**

#### `cells.gene_reads_per_class`

```python
cells.gene_reads_per_class()
```

Total gene reads for each class, weighted by class membership.

This is the soft version of "add up the reads of every cell in a class".
A cell only belongs to a class with some probability $w_{ck}$, so its reads
count in proportion to that probability rather than all-or-nothing:

$$
r_{gk} = \sum_{c=1}^{C} x_{cg}\, w_{ck}
$$

where $x_{cg}$ is the reads of gene $g$ in cell $c$ and $w_{ck}$ is the
probability that cell $c$ is class $k$. This is the weighted total that
mean_gene_reads_per_class then divides by the class size to get an average.

**Returns**

- **`np.ndarray`**
  Shape (G, K): G genes by K cell classes.

#### `cells.mean_gene_reads_per_class`

```python
cells.mean_gene_reads_per_class()
```

Average gene reads for each cell class, in the soft-clustering sense.

Each cell belongs to several classes at once, with probabilities $w_{ck}$,
so the average reads of gene $g$ in class $k$ is a weighted mean over all
cells rather than a plain average:

$$
\overline{r}_{gk} = \frac{\sum_{c=1}^{C} x_{cg}\, w_{ck}}{\sum_{c=1}^{C} w_{ck}}
$$

Here $x_{cg}$ is the reads of gene $g$ in cell $c$ and $w_{ck}$ is the
probability that cell $c$ is class $k$. The top is the weighted read total
for the gene in that class and the bottom is the total probability mass
sitting in the class, so cells count in proportion to how strongly they
belong.

**Returns**

- **`np.ndarray`**
  Shape (G, K): G genes by K cell classes.
