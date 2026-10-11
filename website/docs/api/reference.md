# API reference

::: warning Auto-generated
This page is generated from the pciSeq source by `website/gen_api.py`.
Edit the docstrings in the source, not this file.
:::

Everything here is reachable as `pciSeq.<name>` (plus `VarBayes`, the
model object that [`fit`](#fit) and [`cell_type`](#cell-type) build and
return). The main entry point is [`fit`](#fit). For what the output
DataFrames hold and worked examples of the model attributes, see
[Working with results](./working-with-results).

## `fit`

`pciSeq.app.fit`

```python
fit(*args, **kwargs) -> Tuple[pd.DataFrame, pd.DataFrame]
```

Main entry point for pciSeq cell typing analysis.

``spots`` and ``coo`` may be given as the first two positional arguments or as
keywords. Keyword form is preferred.

**Parameters**

- **`spots`** *(pd.DataFrame)*
  The spots to assign. Needs the columns 'gene_name', 'x' and 'y', plus 'z_plane' for 3D data. Optional 'score' and 'intensity' columns from the spot caller are carried through to geneData; they default to 1.0 and nothing in the model reads them.
- **`coo`** *(list of scipy.sparse.coo_matrix, or numpy array)*
  The label image: each pixel holds the label of the cell it belongs to, and background pixels are labelled 0. Pass it as one sparse matrix per z-plane, or as a dense numpy array of shape (height, width) or (planes, height, width), which is converted to sparse internally. A single coo_matrix is also accepted. More than one plane is treated as 3D. Sparse matrices are modified in place: labels that are not sequential are renumbered, and with `remove_flat_cells` the cells on a single plane are zeroed; pass a copy to keep the original. A dense array is left untouched.
- **`scRNAseq`** *(pd.DataFrame)*
  Cell type definitions: mean expression per gene and cell type, genes as rows and cell types as columns. Required.
- **`opts`** *(dict, optional)*
  Any config values you want to override, e.g. {'max_iter': 500}. See the configuration page for the full list of keys and their defaults.

**Returns**

- **`cellData`** *(pd.DataFrame)*
  One row per cell: the class probabilities, the gene counts and the cell geometry.
- **`geneData`** *(pd.DataFrame)*
  One row per spot: which cell it was assigned to and with what probability.

**Raises**

- **`ValueError`**
  If spots, coo or scRNAseq are missing or malformed, or an option value is out of range.
- **`KeyError`**
  If opts has a key that is not a config option.
- **`TypeError`**
  If an input or an option value has the wrong type.
- **`Exception`**
  Any error raised inside the model is logged and re-raised unchanged. Not converging is not an error: the loop runs to max_iter, logs the convergence status and returns its results.

**Notes**

`spots` and `coo` are given either both as keywords or as the first two positional
arguments. `scRNAseq` and `opts` are keyword only.


## `cell_type`

`pciSeq.app.cell_type`

```python
cell_type(cells: pd.DataFrame, spots: pd.DataFrame, scRNAseq: pd.DataFrame, config: Dict[str, Any]) -> Tuple[pd.DataFrame, pd.DataFrame, VarBayes]
```

Perform cell typing using Variational Bayes algorithm.

**Parameters**

- **`cells`** *(pd.DataFrame)*
  One row per cell, with columns 'label', 'area' and the centroid 'x0', 'y0', 'z0'. stage_data produces it.
- **`spots`** *(pd.DataFrame)*
  One row per spot, with columns 'x', 'y', 'z', 'plane_id', 'gene_name', 'score', 'intensity' and 'label', the label of the cell the spot lies in, 0 for none. stage_data produces it.
- **`scRNAseq`** *(pd.DataFrame)*
  Cell type definitions: mean expression per gene and cell type, genes as rows and cell types as columns. Required.
- **`config`** *(Dict[str, Any])*
  Configuration dictionary containing algorithm parameters

**Returns**

- **`Tuple[pd.DataFrame, pd.DataFrame, VarBayes]`**
  - cellData: DataFrame containing cell typing results - geneData: DataFrame containing gene assignment results - varBayes: The fitted VarBayes model instance

**Raises**

- **`ValueError`**
  If input data is invalid or incompatible
- **`Exception`**
  Any other error raised inside the model is logged and re-raised unchanged. Not converging is not an error, it only logs a warning.


## `stage_data`

`pciSeq.src.preprocess.main.stage_data`

```python
stage_data(spots: pd.DataFrame, coo: List[coo_matrix], cfg: Dict) -> Tuple[pd.DataFrame, Future[Tuple[pd.DataFrame, List[pd.DataFrame]]], pd.DataFrame, Optional[Dict]]
```

Process spots and label images for cell typing analysis.

**Parameters**

- **`spots`** *(pd.DataFrame)*
  Spot data with columns 'gene_name', 'x', 'y', 'z_plane', 'score' and 'intensity'. `fit` fills in 'z_plane', 'score' and 'intensity' when they are missing; a direct call has to supply them.
- **`coo`** *(List[coo_matrix])*
  List of sparse matrices containing cell segmentation. Modified in place: with `remove_flat_cells` the cells that span a single plane are zeroed, and labels that are not sequential are renumbered.
- **`cfg`** *(Dict)*
  Configuration dictionary. Reads 'is3D', 'remove_flat_cells' and 'voxel_size'; writes 'label_map' and 'img_dim'.

**Returns**

- **`cells`** *(pd.DataFrame)*
  Cell properties including position and size
- **`borders_future`** *(Future)*
  Future resolving to (cell_boundaries, cell_boundaries_list). Border extraction runs in the background and only blocks when .result() is called.
- **`processed_spots`** *(pd.DataFrame)*
  Processed spots with cell assignments
- **`label_map`** *(Optional[Dict])*
  Label remapping if labels were reordered


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
application that has its own logging setup, do NOT call this function,
the parent application's handlers will be wiped out. Only call setup_logger()
when pciSeq is the top-level application.

Args:
    level: logging level (e.g. logging.DEBUG, logging.INFO). Defaults to INFO.

Returns:
    The configured 'pciSeq' logger instance.


## `to_spatialdata`

`pciSeq.src.core.utils.spatialdata_export.to_spatialdata`

```python
to_spatialdata(cellData: pd.DataFrame, geneData: pd.DataFrame, coo: List[coo_matrix], varBayes, cfg: Dict)
```

Assemble the run into an in-memory SpatialData object.

Args:
    cellData: cell typing results, one row per cell, original labels.
    geneData: spot results, one row per spot, original labels.
    coo: the segmentation, one sparse plane per z. Either the list fit()
        renumbered in place or a copy with the original labels; which of the
        two it is gets detected, and the store always ends up with the
        original labels.
    varBayes: the fitted model, read for the arrays the two frames do not
        carry (class posterior, spot probabilities, gene panel, reference).
    cfg: the resolved config. voxel_size and label_map are used here.


## `write_spatialdata`

`pciSeq.src.core.utils.spatialdata_export.write_spatialdata`

```python
write_spatialdata(cellData, geneData, coo, varBayes, cfg, out_dir: str) -> str
```

Write the run as a SpatialData zarr store and return its path.


## `add_image`

`pciSeq.src.core.utils.spatialdata_export.add_image`

```python
add_image(store_path: str, img, name: str='background', voxel_size=None, scale_factors=(2, 2, 2)) -> None
```

Add a background image to an existing SpatialData store.

Works like inserting into a database: the store already holds the run
(spots, labels, tables) and this writes one more element into it, without
touching anything else. Kept separate from write_spatialdata because fit()
never sees the image, the same reason stage_image is its own entry point.

Args:
    store_path: path to an existing spatialdata.zarr written by pciSeq.
    img: the image to add. A 2D array (H, W), a 3D stack (Z, H, W), or a
        3D stack with channels (Z, H, W, C). Same shapes stage_image takes.
    name: element name inside the store.
    voxel_size: [x, y, z]. Left as None it is read from the store's own
        provenance, so the microns transform automatically matches the
        other elements. Pass it only for a store that lacks the metadata.
    scale_factors: the multiscale pyramid, each level relative to the one
        before. The default (2, 2, 2) gives four scales. None writes a
        single scale.


## `add_boundaries`

`pciSeq.src.core.utils.spatialdata_export.add_boundaries`

```python
add_boundaries(store_path: str, labels_name: str='cell_labels', name: str='cell_boundaries', voxel_size=None) -> None
```

Add per-plane cell boundary polygons to an existing SpatialData store.

Same insert pattern as add_image, but this one needs no data at all: the
boundaries are derived from the segmentation, and the store already holds
the segmentation as the labels element. They are extracted here with the
same chain code tracing the pipeline uses, so what goes in matches what
pciSeq would have drawn.

Shapes in SpatialData are strictly 2D, so a 3D run gets one shapes element
per plane, cell_boundaries_plane_000 and so on, each a set of polygons
indexed by the cell label. A cell spanning 12 planes appears as 12
polygons, one per element, all under its own label, which is exactly how
pciSeq thinks of boundaries anyway. A 2D run gets a single element.

Args:
    store_path: path to an existing spatialdata.zarr written by pciSeq.
    labels_name: the labels element to trace.
    name: element name, used as a prefix on 3D runs.
    voxel_size: [x, y, z]. Left as None it is read from the store's own
        provenance, same as add_image.


## `read_tiles`

`pciSeq.src.tiling.read_tiles.read_tiles`

```python
read_tiles(mbtiles: str, plane: int=0, bbox: Optional[Tuple[float, float, float, float]]=None, width: Optional[int]=None, zoom: Optional[int]=None) -> Tuple[Image.Image, float]
```

Read part of a plane out of an mbtiles pyramid.

**Parameters**

- **`mbtiles`** *(str)*
  Path to the .mbtiles file, as written by :func:`pciSeq.stage_image`.
- **`plane`** *(int, default 0)*
  The z-plane to read. For a 2D image there is only plane 0.
- **`bbox`** *(tuple of float, optional)*
  The region to read, ``(x0, y0, x1, y1)`` in the pixel coordinates of the original image, the same coordinates as `cellData` and the spots. The box is clamped to the image. The default is the whole plane.
- **`width`** *(int, optional)*
  Width of the returned image in pixels. The pyramid level is picked to cover it and the result is resampled to exactly this width. The default returns the region at its original scale, one pixel per image pixel.
- **`zoom`** *(int, optional)*
  Read this pyramid level instead of choosing one. Mostly useful for inspecting the file itself.

**Returns**

- **`image`** *(PIL.Image.Image)*
  The region, in RGB.
- **`scale`** *(float)*
  Zoom factor, returned width over box width. A point ``(x, y)`` of the original image is at ``((x - x0) * scale, (y - y0) * scale)`` in the returned one.

**Notes**

This is not a lossless recovery of the original data. The tiles are JPEG and every
level was produced by resizing, so the crop is a close visual copy, not the raw
pixels of the image the pyramid was built from. Use it for figures and for checking
a segmentation against the image, not for measurements.

`width` sets the size of the result and therefore how much detail is read. The
pyramid level is chosen as the cheapest one that covers the request, and the result
is resampled to exactly `width`. Levels only go as high as the pyramid does, so a
`width` beyond the top level is allowed and enlarges the result, which is
magnification, not detail.

**Examples**

A cell and its surroundings, one pixel per image pixel:

>>> im, scale = read_tiles('dapi.mbtiles', plane=57, bbox=(5393, 702, 5603, 842))
>>> im.size, scale
((210, 140), 1.0)

The same region as a 1200 pixel wide panel. The box is 210 across, so `scale` is
5.71 and the panel is magnified: the level read is the top of the pyramid, which
renders that box 10.2 times up, and it is resampled down to 5.71.

>>> im, scale = read_tiles('dapi.mbtiles', plane=57, bbox=(5393, 702, 5603, 842),
...                        width=1200)

A whole 6408 by 4382 plane, trimmed to 3:2 and shown 1200 wide. Here `scale` is
1200 / 6408 = 0.187, and a cell at (5498.2, 772.5) lands at
``((5498.2 - 0) * 0.187, (772.5 - 55) * 0.187)``, so about (1029, 134) in the
returned image. `bbox=None` would give the untrimmed plane.

>>> im, scale = read_tiles('dapi.mbtiles', plane=57, bbox=(0, 55, 6408, 4327),
...                        width=1200)


## `stage_image`

`pciSeq.src.tiling.stage_image.stage_image`

```python
stage_image(img, out_dir=None, zoom_levels=None, name=None, description=None, plane_prefix='plane_', use_buffer=True, tint=None, progress=True)
```

Turn an image (or z-stack) into an MBTiles file the viewer can read.

This builds the tiled, multi-resolution background that pciSeq Viewer uses
as its slippy-map base layer. Give it a single 2D image or a whole 3D stack
and it writes one `.mbtiles` file.

**Parameters**

- **`img`** *(np.ndarray or str)*
  The image to tile: a 2D array (H, W), a 3D stack (Z, H, W), a 3D stack with channels (Z, H, W, C), or a path to a 2D image file (legacy).
- **`out_dir`** *(str, optional)*
  Directory for the `.mbtiles` file. Defaults to the system temp directory.
- **`zoom_levels`** *(int, optional)*
  The deepest zoom level. Levels 0 to `zoom_levels` are written, and level z is 256 * 2**z pixels on its longer side. If not given, it is the first level that holds the image at full resolution, for example 5 (8192 pixels) for an image 6408 pixels wide. A deeper level adds no detail, only file size; the viewer enlarges the deepest level when zooming in past it.
- **`name`** *(str, optional)*
  Short identifier for the dataset. Also used as the output filename, e.g. `name="S10_gcamp_10"` writes `S10_gcamp_10.mbtiles`. If empty, the file is named `output.mbtiles`.
- **`description`** *(str, optional)*
  Longer description of the dataset.
- **`plane_prefix`** *(str, optional)*
  Prefix for the per-plane names. Default is "plane_".
- **`use_buffer`** *(bool, optional)*
  If True (the default) the tiles are built in memory and inserted straight into the MBTiles database. If False they are written to disk first, which uses less memory but more disk I/O.
- **`tint`** *(str, optional)*
  Hex colour like "#00FF00" the viewer uses to tint this grayscale layer. If omitted, the layer is shown in plain grayscale.
- **`progress`** *(bool, optional)*
  If True (the default) show two per-plane tqdm bars, one for tiling and one for the db writing. Set to False for headless/quiet runs.

**Returns**

- **`str`**
  Path to the created `.mbtiles` file.

**Notes**

Requires libvips. If it is not installed, `pciSeq.stage_image()` falls back
to a stub that only logs a warning.


## `tile_maker`

`pciSeq.src.tiling.stage_image.tile_maker`

```python
tile_maker(img, zoom_levels=None, out_dir='./tiles', plane_prefix='plane_', progress_bar=None)
```

Makes a pyramid of tiles from an image.

Args:
    img: One of:
        - str: path to a 2D image file (TIFF, PNG, JPEG, etc.)
        - numpy array (H, W): single 2D grayscale image
        - numpy array (Z, H, W): 3D stack of grayscale images
        - numpy array (Z, H, W, C): 3D stack with channels
    zoom_levels: (int) The deepest zoom level. Levels 0 to zoom_levels are written. If not
        given it is the first level big enough to hold the image, see stage_image.
    out_dir: (str) Output folder for the tile pyramid. Will be deleted and recreated if exists.
    plane_prefix: (str) Prefix for plane subdirectories when processing 3D images.
                  Default is "plane_" resulting in "plane_0", "plane_1", etc.
    progress_bar: (tqdm, optional) if given, ticked once per plane instead of
                  logging a per-plane line. Used by stage_image to drive its bars.

Returns:
    dict with keys:
        - 'original_dims': [width, height] of the original input image
        - 'num_planes': number of planes processed
        - 'zoom_levels': the deepest zoom level that was written


## `open_run`

`pciSeq.src.mcp.tools.open_run`

```python
open_run(path)
```

Open a finished run. The entry point every other tool goes through.


## `VarBayes`

`pciSeq.src.core.main.VarBayes`

```python
VarBayes(cells_df: pd.DataFrame, spots_df: pd.DataFrame, scRNAseq: pd.DataFrame, config: Dict[str, Any])
```

The variational Bayes model: assigns spots to cells and cells to classes.

**Parameters**

- **`cells_df`** *(pd.DataFrame)*
  One row per cell, with columns 'label', 'area' and the centroid 'x0', 'y0', 'z0'. stage_data produces it.
- **`spots_df`** *(pd.DataFrame)*
  One row per spot, with columns 'x', 'y', 'z', 'plane_id', 'gene_name', 'score', 'intensity' and 'label', the label of the cell the spot lies in, 0 for none. stage_data produces it.
- **`scRNAseq`** *(pd.DataFrame)*
  Cell type definitions: mean expression per gene and cell type, genes as rows and cell types as columns.
- **`config`** *(dict)*
  The resolved configuration.

**Attributes**

- <a id="cells-spots-genes-single-cell-celltypes"></a>**`cells, spots, genes, single_cell, cellTypes`** *(object)*
  The parts of the model. Working with results describes the arrays on each.
- <a id="nc-ns-ng-nk"></a>**`nC, nS, nG, nK`** *(int)*
  Number of cells (including the background row 0), spots, genes and classes (including Zero).
- <a id="config"></a>**`config`** *(dict)*
  The configuration of the run: the defaults, the `opts` overrides, and the runtime keys `is3D`, `img_dim` and `label_map`.
- <a id="metadata"></a>**`metadata`** *(dict)*
  Provenance recorded when the model is built: `version`, `branch`, `commit`, `commit_date`, `created_at`, `python_version`, `os` and `package_versions`.
- <a id="has-converged"></a>**`has_converged`** *(bool)*
  True when the loop stopped because the change fell below `CellCallTolerance`, False when it ran to `max_iter`.
- <a id="iter-delta"></a>**`iter_delta`** *(list of float)*
  The largest change in the spot assignment probabilities, one entry per iteration.
- <a id="iter-num"></a>**`iter_num`** *(int)*
  Index of the last iteration run, counting from 0.

::: tip Obtaining a fitted instance
`VarBayes` is not instantiated directly in normal use. [`fit`](#fit) and [`cell_type`](#cell-type) construct and run it. `cell_type` returns the fitted instance; `fit` does not, but when `save_data=True` (the default) the fitted model is serialised to `<output_path>/pciSeq/data/debug/pciSeq.pickle` (`output_path` defaults to a temporary directory). The attributes and methods below operate on a loaded instance; [Working with results](./working-with-results) walks through the main ones with examples.

```python
import pandas as pd

obj = pd.read_pickle('<output_path>/pciSeq/data/debug/pciSeq.pickle')

obj.metadata
obj.check_cell(my_label=42, user_class='Astro')
```
:::

### Methods

#### `check_spot`

```python
check_spot(spot_id, show_plot=True)
```

Break down the assignment of a spot to its candidate cells.

The score of each candidate cell is the sum of the terms from the last
spot-to-cell update: the spatial log-likelihood ('Gaussian fit'), the expected log mean
expression under the cell's class probabilities ('class expression'), the expected
log theta ('cell scale'), the expected log gamma ('cell-gene scale'), the log eta
of the spot's gene ('gene efficiency') and the inside-cell bonus ('bonus'). The
score of the background is the log misread density of the gene. Probabilities are
the softmax of the scores.

**Parameters**

- **`spot_id`** *(int)*
  Spot id, the index of the spots table.
- **`show_plot`** *(bool, default True)*
  Draw the score and probability charts.

**Returns**

- **`pd.DataFrame`**
  One row per candidate cell, indexed by cell label, and a final 'background' row. Columns are the score terms above, 'misread', 'sum' and 'prob'.

#### `check_cell`

```python
check_cell(my_label, user_class, top_n=10, show_plot=True, top_classes=5)
```

Compare the assigned cell type of a cell with another cell type.

The assigned type is the type with the highest probability in classProb. The
per-gene log-likelihoods, log prior and spatial (MRF) term are those used in
the last cell type update. The figure shows the genes that most favour each
type, the three score components for both types, and the posterior over all
types for the most likely ones.

**Parameters**

- **`my_label`** *(int)*
  Cell label, as in the input segmentation.
- **`user_class`** *(str)*
  Cell type to compare against. Must differ from the assigned type.
- **`top_n`** *(int, default 10)*
  Maximum number of genes shown on each side: the genes that most favour the assigned type and those that most favour user_class. Only genes with a nonzero difference in that direction are shown, so a side can have fewer.
- **`show_plot`** *(bool, default True)*
  Draw the figure.
- **`top_classes`** *(int, default 5)*
  Number of types shown in the posterior chart, ordered by probability. user_class is added if it is not among them.

**Returns**

- **`gene_expression_data`** *(pd.DataFrame)*
  One row per selected gene. Columns are the mean counts of the gene in cells of each type, the expected count of the gene in this cell under each type, and its counts in this cell. The counts are sums of spot assignment probabilities, so they are estimates, not measurements.
- **`contr`** *(pd.DataFrame)*
  Log-likelihood of each gene in this cell under the two types, and their difference, 'diff'.
- **`fig`** *(matplotlib.figure.Figure or None)*
  The figure, or None if show_plot is False.

**Raises**

- **`ValueError`**
  If user_class is not a cell type, or is the assigned type of the cell.

#### `to_internal`

```python
to_internal(label)
```

Convert a segmentation cell label to the internal one.

The internal label is the cell's row index in the arrays in the pickle file,
1 to nC - 1, with row 0 the background. It differs from the segmentation label only when the
input labels were not sequential, which usually happens because `remove_flat_cells`
dropped the cells that span a single plane.

**Parameters**

- **`label`** *(int or list of int)*
  Cell label, as in the input segmentation.

**Returns**

- **`int or list of int`**
  The row to index the arrays with. The same value back when no renumbering took place.

**Raises**

- **`KeyError`**
  If the label is not in the segmentation.

#### `to_external`

```python
to_external(label)
```

Convert an internal cell label to the segmentation one.

The opposite of `to_internal`. Use it on a row index taken from an array in the
pickle file or from diagnostics.db, to get the label the segmentation gave that cell.

**Parameters**

- **`label`** *(int or list of int)*
  Row index into the arrays in the pickle file.

**Returns**

- **`int or list of int`**
  The segmentation label. The same value back when no renumbering took place.

**Raises**

- **`KeyError`**
  If the row index is out of range.

#### `read_tsv`

```python
read_tsv(filepath)
```

Read a tsv file written by pciSeq into a DataFrame.

Columns that hold lists or dicts are parsed back from text.

**Parameters**

- **`filepath`** *(str)*
  Path to the file, e.g. cellData.tsv or geneData.tsv.

**Returns**

- **`pd.DataFrame`**

#### `heatmap_counts_per_class`

```python
heatmap_counts_per_class()
```

Display the interactive heatmap.

#### `cells.gene_reads_per_class`

```python
cells.gene_reads_per_class()
```

Calculate total (weighted by class prob) gene reads for each class.

Returns:
    np.ndarray: Shape (G, K) total reads per class and gene

#### `cells.mean_gene_reads_per_class`

```python
cells.mean_gene_reads_per_class()
```

Calculate the average gene reads for each cell class/type in a soft clustering setup.

In soft clustering, each cell belongs to multiple classes with probabilities $w_{ck}$.
The average number of reads for gene $g$ in class $k$ is computed as:

$$\overline{r}_{gk} = \frac{\sum_{c=1}^{C} x_{cg} \cdot w_{ck}}{\sum_{c=1}^{C} w_{ck}}$$

Where:
    - $x_{cg}$: Number of reads for gene $g$ in cell $c$
    - $w_{ck}$: Probability that cell $c$ belongs to class $k$
    - The numerator is the total weighted sum of reads for gene $g$ in class $k$
    - The denominator is the total probability mass of class $k$

Returns:
    np.ndarray: Shape (G, K), where:
        G = number of genes
        K = number of cell classes/types


## `Run`

`pciSeq.src.mcp.tools.Run`

```python
Run(path)
```

A finished pciSeq run, opened for questions.

**Parameters**

- **`path`** *(str or Path)*
  The run folder, or the diagnostics.db itself. Anything above the database works, the file is searched for underneath.

**Notes**

diagnostics.db keeps `gamma_assigned`, the gamma of the winning class only, not
the whole (nC, nG, nK) array. Questions about gamma under a different class need
the pickle.

::: tip Obtaining a Run
[`open_run`](#open-run) builds one from the output folder. Many of the methods below are also functions the agent can call, under the same names. The [MCP server](./mcp-server) page lists the agent's functions.

```python
from pciSeq.src.mcp.tools import open_run

run = open_run('<output_path>')
run.explain_cell(2413)
```
:::

### Methods

#### `summary`

```python
summary()
```

What this run is: its size, the code that made it, and what it can answer.

#### `cell`

```python
cell(label)
```

The headline facts about one cell.

#### `explain_cell`

```python
explain_cell(label, vs_class=None, top_n=10)
```

Why this cell got its class, gene by gene.

The class score is the sum of the per gene negative binomial contributions,
plus the class prior, plus the spatial term. This recomputes those from the
stored arrays and compares the assigned class against another one.

The recomputation uses the final eta, and the loop updates eta after the last
class update, so it reproduces the stored probabilities exactly only once the
run has converged. On an unconverged run the argmax still agrees but the
probabilities can be off in the second decimal.

**Parameters**

- **`label`** *(int)*
  Cell label, as in your segmentation.
- **`vs_class`** *(str, optional)*
  The class to compare against. Defaults to the runner up.
- **`top_n`** *(int, default 10)*
  How many genes to report on each side.

#### `explain_spot`

```python
explain_spot(spot_id)
```

Why this spot went to the cell it did, term by term.

One row per candidate cell plus the background. The terms are the ones the
model used in its last spot update, so the probabilities here are the ones
in geneData, at full float32 precision rather than the 3 decimals the file
keeps. Use this rather than spot_row when a probability below 0.0005 matters.

#### `cell_counts`

```python
cell_counts(label, gene=None)
```

How many reads a cell holds. Soft, and the answer says so.

#### `spots_in_cell`

```python
spots_in_cell(label, gene=None)
```

How many spots physically sit inside a cell's segmentation mask.

Hard containment, no probabilities. A spot can sit outside every cell and
still be assigned to one, so this is a different number from `cell_counts`.

#### `spots_of_cell`

```python
spots_of_cell(label, min_prob=None, gene=None)
```

Which spots belong to a cell, under one of two definitions.

With `min_prob` unset: the spots whose most likely parent is this cell, the
argmax. With `min_prob` set: every spot with a probability on this cell above
it, which is what cellData.spot_id holds at 0.0001. The two are different
lists, and the answer says which one it is. Either way every spot comes with
its probability, sorted highest first. `gene` keeps only that gene's spots,
so the count of them is the tool's number and not the agent's.

#### `cell_row`

```python
cell_row(label)
```

The cellData row of one cell, value for value.

Everything the tsv has except the two drawing columns. Read straight out of
cellData.tsv when it is there, which is both exact and fast, a targeted grep
on a 250 MB file takes about 30 ms. Without the tsv the row is rebuilt from
the viewer files, which matches on every column except the boundary cases of
spot_id; see `_exact_spot_ids`.

#### `spot_row`

```python
spot_row(spot_id)
```

The geneData row of one spot, value for value.

Read out of geneData.tsv when it is there, otherwise rebuilt from the viewer
files; `source` says which. Both give the same numbers, the viewer files are
written from the same frame.

neighbour_prob is rounded to 3 decimals here because that is how geneData
stores it (summary.py:102), so a candidate holding 0.0004 shows up as 0.0 and
you cannot tell it apart from one holding nothing at all. That is fine for
reporting what the file says, and it is why cell_row's spot_id list can
contain spots whose probability reads 0.0: the cut-off there is 0.0001,
applied before the rounding. When the small probabilities matter use
explain_spot, which recomputes them from diagnostics.db at full precision.

#### `cell_image`

```python
cell_image(label, context=False, plane=None, width=1200, channel=None, neighbours=False, mbtiles=None)
```

A picture of one cell on the background image of the run.

The background is stitched from the viewer's tile pyramid with read_tiles, so
it is a close copy of the image, not the raw pixels (the tiles are jpeg). Two
kinds of picture:

* close-up (context=False): the cell outlined in red and every other cell on
  that plane in blue, in a 3:2 window about four and a half times the size of
  the cell. With neighbours=True only the cells the mrf term listens to are
  outlined, the picture for 'why did its neighbours make it this class'. Those
  on another plane have no outline here and the answer lists them.
* context (context=True): the whole plane trimmed to 3:2 with a ring round the
  cell, to show where in the tissue it sits.

The plane is the one given, else the plane of the cell centroid when the run
carries voxel_size, else the plane where the cell outline is biggest. Returns
the image (PIL) and a dict of what was drawn. The MCP tool also takes
save_as; here just call .save on the image.

channel picks the background image when the run has more than one (DAPI,
GCaMP, ...); with several and no channel it refuses and lists them. With one
image channel is not checked, the one image is used.

#### `plane_image`

```python
plane_image(plane=None, bbox=None, width=1200, channel=None, mbtiles=None)
```

The background image of one plane, whole or a part of it, with nothing
drawn on top.

Stitched from the viewer's tile pyramid with read_tiles, so like cell_image
it is a close copy of the image, not the raw pixels. The plane defaults to
the middle one of the stack. bbox is (x0, y0, x1, y1) in image pixels, the
same coordinates as the cells and spots, and is clamped to the image; leave
it out for the whole plane, untrimmed. Returns the image (PIL) and a dict of
what was read. The MCP tool also takes save_as; here just call .save on the
image. channel works as in cell_image.

#### `to_internal`

```python
to_internal(label)
```

Segmentation label to the row it occupies in the arrays.

#### `to_external`

```python
to_external(row)
```

The row of an array back to the segmentation label.
