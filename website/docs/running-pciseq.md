# Running pciSeq

pciSeq is run from Python by calling [`fit`](api/reference.md#fit) with the spots, the
segmentation and the cell type definitions, or from the
[command line](api/command-line.md) with the same inputs listed in a config file.

## From Python

```python
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import pciSeq

spots = pd.read_csv("spots.csv")                          # gene_name, x, y, z_plane
coo = [coo_matrix(plane) for plane in np.load("masks.npy")]  # one label matrix per plane
scRNAseq = pd.read_csv("scRNAseq.csv", index_col=0)        # genes x cell types

opts = {
    "voxel_size": [0.28, 0.28, 0.7],
    "exclude_genes": ["Zbtb20", "Ddit4l"],
    "output_path": "out/run1",
}

cellData, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=scRNAseq, opts=opts)
```

`fit` returns two DataFrames, one row per cell and one row per spot, described in
[Working with results](api/working-with-results.md). With `save_data` on, the default,
the same results are also written under `<output_path>/pciSeq/data/`.

## Inputs

| argument | type | required | content |
| --- | --- | --- | --- |
| `spots` | DataFrame | yes | One row per spot. Columns `gene_name`, `x`, `y`, and `z_plane` for 3D data. `score` and `intensity` are optional; they are carried through to the output and default to 1.0. |
| `coo` | list of `scipy.sparse.coo_matrix`, or a numpy array | yes | The segmentation as a label image: each pixel holds the label of the cell it belongs to, and background pixels are labelled 0. Pass it as one sparse matrix per z-plane, or as a dense numpy array of shape `(h, w)` or `(planes, h, w)`, which is converted to sparse internally. |
| `scRNAseq` | DataFrame | yes | Mean expression per gene and cell type, genes as rows and cell types as columns, gene names in the index. Genes absent from the panel are dropped from the spots with a warning. |
| `opts` | dict | no | Settings to override, see [Configuration](api/configuration.md). Anything left out keeps its default. |

Coordinates in `spots` are in pixels of the label image. Spots and segmentation have to
be in the same frame.

## 2D and 3D

The dimensionality is read off `coo`. A single matrix, or a list with one element, is 2D.
A list with more than one plane is 3D.

- **2D.** `z_plane` may be left out of `spots`, it is set to 0. The legacy column name
  `Gene` is accepted for `gene_name`. `voxel_size` stays at its default `[1, 1, 1]`.
- **3D.** `z_plane` is required and indexes the plane in `coo`. `voxel_size` should be
  set to the physical size of a voxel as `[x, y, z]`, since the z step is usually
  coarser than the pixel and the distances along z are wrong without it.
  `remove_flat_cells`, on by default, drops cells that appear on a single plane.

## Settings

Every setting and its default is in `pciSeq/config.py`; the
[configuration page](api/configuration.md) is generated from it. `opts` overrides
individual keys, the rest keep their defaults.

The settings are of two kinds:

- Experiment settings describe the data and follow from it: `voxel_size`,
  `exclude_genes`, `output_path`, the viewer options.
- Hyperparameters set the behaviour of the model. The user picks them, and
  usually tries a few values before settling on one. The next section goes
  through them.

## Fine-tuning

The hyperparameters below are listed in the order in which they are typically adjusted.

| Setting | Controls |
| --- | --- |
| [`Inefficiency`](api/configuration.md#inefficiency) | The ratio of the sensitivities of the two technologies: the fraction of a cell's transcripts the spatial assay detects, over the fraction scRNA-seq detects. It scales the expected counts of every class. At 0.2, for example, the spatial assay is taken to detect one fifth of what scRNA-seq reports. |
| [`rTheta`](api/configuration.md#rtheta) | Each cell has a scale factor $\theta_c$ that rescales the expected counts of every class to the cell's own total count; `rTheta` sets how strongly $\theta_c$ is held at 1. A low value, towards 1, lets the data set $\theta_c$; a high value, in practice an order of magnitude above the typical number of counts in a cell, pushes $\theta_c$ to 1, so that it no longer changes the expected counts. A low value lets a near-empty cell take a real class on one or two spots; see [the prior strength](the-model/scale-factors.md#r-theta). |
| [`mrf_beta`](api/configuration.md#mrf-beta) | Spatial regularization strength. Larger values pull a cell's class towards the classes of its neighbours. The bonus a class receives is `mrf_beta` times the sum of the neighbours' probabilities for that class, each neighbour weighted by its distance; the weights sum to `nNeighbors`. `0` disables the term. Useful where cell types are strongly localised, in layers or regions, as in the hippocampus. |
| [`MisreadDensity`](api/configuration.md#misreaddensity), [`rRho`](api/configuration.md#rrho) | The prior mean of the per-gene background density, in spots per unit volume, and the strength of that prior. `rRho` counts as that many spots against the gene's background count, so at the default of 1 the prior has very little effect and the density is set by the data. |
| [`cell_type_weights`](api/configuration.md#cell-type-weights) | The prior probability of each class. Types not listed share the remaining probability equally, so the default `{'Zero': 0.5}` gives Zero half the prior mass and splits the other half among the real types. There is no `default` key. |
