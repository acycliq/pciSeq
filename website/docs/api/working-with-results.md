---
title: Working with results
description: What pciSeq.fit() returns, the columns of cellData and geneData, and how to inspect the fitted model.
---

# Working with results

The two DataFrames [`fit`](./reference#fit) returns, their columns, and the fitted
model object behind them.

## Return values

```python
cellData, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=ref)
```

`cellData` has one row per cell; `geneData` has one row per spot. Several columns
hold **parallel lists**: two columns of the same length whose entries line up by
position, so element `i` of one matches element `i` of the other.

The notebook
[`1_pciSeq.ipynb`](https://github.com/acycliq/pciSeq_3d/blob/dev_3d/notebooks/1_pciSeq.ipynb)
runs `fit` on a 2D section of mouse hippocampus and goes through both DataFrames.

## cellData

One row per cell. The background pseudo-cell is dropped, so every row is a real
cell.

| Column | Type | Meaning |
| --- | --- | --- |
| `Cell_Num` | int | Cell label. The cell's identity, matching the label image. |
| `X`, `Y`, `Z` | float | Fitted cell centroid, the centre of the cell's Gaussian after the [cell geometry](../the-model/gaussian.md) update. `Z` is 0 for 2D data. |
| `X0`, `Y0`, `Z0` | float | The cell's centroid as given by the segmentation. It can differ from `X`, `Y`, `Z`; how much depends on the prior on the cell's centre, [`cell_centroid_prior`](./configuration.md#cell-centroid-prior), and on the number of spots in the cell. `Z0` is 0 for 2D data. |
| `ClassName` | list of str | Candidate cell classes, sorted by descending probability. |
| `Prob` | list of float | Probability of each class, lined up with `ClassName`. |
| `Genenames` | list of str | Genes assigned to the cell, sorted by descending count. |
| `CellGeneCount` | list of float | Expected count of each gene, lined up with `Genenames`: the sum of the assignment probabilities of the gene's spots, see [probabilistic output](../index.md#probabilistic-output). |
| `spot_id` | list of list of int | The ids of the spots that make up each gene's count, lined up with `Genenames`: the gene's spots whose probability of belonging to this cell is above 0.0001. Their probabilities add up to the count. |
| `gaussian_contour` | list | The 3-sigma ellipse outline of the cell, for drawing. Not written to `cellData.tsv`. |
| `sphere_scale`, `sphere_rotation` | list | The fitted covariance as an ellipsoid: the three semi-axes of the 3-sigma ellipsoid, longest first, and its Euler angles in radians. `ellipsoid_to_cov` rebuilds the covariance from them. 3D data only. |

::: tip Coordinate units and the z-plane
`X`, `Y` and `Z` are isotropic coordinates in units of the x-voxel
(`voxel_size[0]`), not microns. pciSeq rescales every axis to the x-voxel size so
Euclidean distances are comparable in all directions. With the default isotropic
`voxel_size` of `[1, 1, 1]` the values are plain pixels and `Z` equals the
z-plane index. With anisotropic voxels, `Z` is the plane index scaled by
`voxel_size[2] / voxel_size[0]`, so to recover the (fractional) plane of a
fitted cell centroid:

```python
plane = Z * voxel_size[0] / voxel_size[2]
```

Spots in `geneData` carry the integer plane directly in `plane_id`, so they need
no conversion.
:::

### cellData column relationships

`ClassName` and `Prob` are a parallel pair. The class typing of a cell is read
off them together: `ClassName[0]` is the most likely class and `Prob[0]` is its
probability, `ClassName[1]` the runner-up, and so on. Only classes with
probability above `0.001` are kept, so the lists are usually short.

```python
row = cellData.iloc[0]
dict(zip(row['ClassName'], row['Prob']))
# {'Astro': 0.71, 'Oligo': 0.22, 'Zero': 0.05, ...}
```

`Genenames`, `CellGeneCount` and `spot_id` are a parallel triple, in the same
way: `Genenames[i]` was seen `CellGeneCount[i]` times in this cell, contributed
by the spots in `spot_id[i]`.

```python
dict(zip(row['Genenames'], row['CellGeneCount']))
# {'Plp1': 14.0, 'Mbp': 9.0, ...}
```

## geneData

One row per spot, recording where each spot was assigned.

| Column | Type | Meaning |
| --- | --- | --- |
| `gene_name` | str | The gene the spot belongs to. |
| `gene_id` | int | Integer index of the gene into the gene panel. |
| `spot_id` | int | The spot's identity. |
| `x`, `y`, `z` | float | Spot coordinates. `z` is 0 for 2D data. |
| `plane_id` | int | The z-plane the spot sits on. |
| `neighbour` | int | The most likely parent cell label (the top of `neighbour_array`). |
| `neighbour_array` | list of int | The candidate parent cells, sorted by descending probability. The background cell `0` means a misread. |
| `neighbour_prob` | list of float | Probability of each candidate, lined up with `neighbour_array`. |
| `inside_cell` | int | The cell whose segmentation mask the spot falls in, as a `Cell_Num`, `0` for none. Physical containment, not an assignment: `neighbour` is the cell the model chose, and a spot can sit outside every cell and still be assigned to one. |
| `omp_score`, `omp_intensity` | float | Spot detection score and intensity. 1.0 when the input spots had no such columns. |
| `is_hard_misread` | uint8 | 1 when the argmax over the candidate probabilities is the background (`neighbour` is 0). In the saved files only, not in the returned DataFrame. |

### geneData column relationships

`neighbour_array` and `neighbour_prob` are a parallel pair: `neighbour_array[i]`
is a candidate cell and `neighbour_prob[i]` is the probability the spot belongs
to it. `neighbour` is `neighbour_array[0]`, the most probable candidate.

The background cell, label `0`, is always one of the candidates. Its probability
is the chance the spot is a misread, so a spot with a high probability on `0` was
most likely noise.

```python
row = geneData.iloc[0]
dict(zip(row['neighbour_array'], row['neighbour_prob']))
# {458: 0.88, 12: 0.09, 0: 0.03}   # cell 458 wins; 0 is the misread chance
```

## Saved files

With `save_data` on (the default), the same results are written under
`<output_path>/pciSeq/data/`:

```
data/
    tsv/
        cellData.tsv
        geneData.tsv
        cellBoundaries.tsv
    viewer_data/
        arrow_spots/          one feather file per plane
        arrow_cells/
        arrow_boundaries/     one feather file per plane
        diagnostics/diagnostics.db
    spatialdata.zarr/
    debug/pciSeq.pickle
```

`tsv/` holds the two DataFrames and the cell outlines. `viewer_data/` is what
[pciSeq Viewer](../viewer/overview.md) reads. `spatialdata.zarr` is described on the
[SpatialData store](./spatialdata-store.md) page, and `debug/pciSeq.pickle` is the fitted model, see below.

## The diagnostics.db file

`viewer_data/diagnostics/diagnostics.db` is an SQLite database. It holds, for every cell
and every spot, the terms from which its assignment was computed, together with the
settings of the run. pciSeq Viewer and the agent read it.

| Table | Rows | Contents |
| --- | --- | --- |
| `metadata` | one per key | The settings (`config`), the convergence record (`run`), the version of pciSeq (`pciSeq_provenance`), the class names, the gene panel and the per-gene quantities. Values are text; lists and dictionaries are stored as JSON. |
| `cells` | one per cell, plus the background | The class probabilities, the scale factor, the spatial term, the neighbouring cells, the centre and covariance, and the gene log-likelihoods. |
| `spots` | one per spot | The candidate cells of the spot and the terms of its score for each. |
| `docs` | one per page | The documentation of the version of pciSeq that made the run. |

The `metadata` table can be read directly:

```python
import json
import sqlite3

con = sqlite3.connect('<output_path>/pciSeq/data/viewer_data/diagnostics/diagnostics.db')
value, = con.execute("select value from metadata where key = 'config'").fetchone()
config = json.loads(value)
config['mrf_beta']
# 1.0
```

Most columns of `cells` and `spots` are arrays stored as bytes, and their rows are
indexed by the [internal label](#cell-identifiers). The layout of these two tables may
change between versions of pciSeq. The same terms are available from the fitted model
through [`check_cell` and `check_spot`](../explaining-the-calls/overview.md), and from
the [agent](../accessing-the-ai-agent.md).

## Inspecting the fitted model

`cellData` and `geneData` are summaries. The full state of the run is on the fitted
`VarBayes` model. [`VarBayes`](./reference#varbayes) lists its attributes and methods;
this section covers the main ones.

[`cell_type`](./reference#cell-type) returns the model directly. `fit` does not,
but with `save_data=True` (the default) it serialises the model to
`<output_path>/pciSeq/data/debug/pciSeq.pickle`:

```python
import pandas as pd
obj = pd.read_pickle('<output_path>/pciSeq/data/debug/pciSeq.pickle')
```

The model holds the run dimensions and a small object graph: `obj.cells`,
`obj.spots`, `obj.genes`, `obj.single_cell`, `obj.cellTypes`.

```python
obj.nC, obj.nS, obj.nG, obj.nK   # cells, spots, genes, classes
```

### Cell identifiers

A cell can be identified by two numbers. The **segmentation label** is the one it carries in the label
image passed to `fit`. The **internal label** is its row index in the arrays in the pickle
file, `1` to `nC - 1`, with row `0` the background.

The two are the same number unless the segmentation labels were not sequential, in which
case pciSeq renumbers them. For 3D data that is the usual case rather than the exception,
since [`remove_flat_cells`](./configuration.md#remove-flat-cells), on by default, drops the
cells that appear on a single plane.

When pciSeq saves the data then the cell identifiers have been mapped back to the original
segmentation label. Hence the rule to remember is:
1. the pciSeq exported flatfiles (tsv, arrow etc) are in segmentation labels, and so is the viewer
2. the arrays in the pickle file, and `diagnostics.db`, are indexed by the internal label;
   in `diagnostics.db` the column is called `internal_label`

::: warning Looking up a cell in the pickle arrays
Suppose we are interested in cell 18223, the label the segmentation uses and the viewer
shows. Its class probabilities are in `cells.classProb`, which is indexed by the internal
label, so 18223 has to be converted first:

```python
obj.cells.classProb[obj.to_internal(18223)]   # correct
obj.cells.classProb[18223]                    # a different cell, or IndexError
```

The methods of the model, such as [`check_cell`](./reference.md#check-cell), take the
segmentation label and convert it themselves, so this applies only to indexing an array
directly.
:::

[`to_internal`](./reference.md#to-internal) and
[`to_external`](./reference.md#to-external) convert either way, and return what they were
given when no renumbering took place.

### Class probabilities (`cells.classProb`)

`cells.classProb` is a `(nC, nK)` array: one **row per cell**, one **column per
class**, and each row sums to 1. Row `0` is the background. The
columns are labelled by `cells.class_names`, the last one `Zero`. As a DataFrame:

```python
probs = pd.DataFrame(obj.cells.classProb, columns=obj.cells.class_names)
probs.iloc[1:].idxmax(axis=1)   # most likely class for each real cell
```

Rows follow the internal cell labels described above.

```python
obj.cells.classProb[obj.to_internal(18223)]
```

`cellData`'s `ClassName` and `Prob` columns are built from this array: for each
cell the classes are sorted by descending probability and any below `0.001` are
dropped. `cells.classProb` is the raw form, keeping every class in the fixed
`class_names` order.

### Spot assignment (`spots.parent_cell_prob`)

`spots.parent_cell_prob` is `(nS, nNeighbors + 1)`: for each spot, the
probability over its nearest candidate cells. The matching cell labels are in
`spots.parent_cell_id`, in the same internal numbering. The **last column is the
background**, i.e. the misread probability. `geneData`'s `neighbour_array` and
`neighbour_prob` are the same numbers sorted by descending probability; these two
arrays are the raw form, with the candidates in a fixed order.

### Genes (`genes.gene_panel`, efficiency, misread density)

`genes.gene_panel` is the master list of genes, length `nG`. Every gene-indexed
array below is in this order.

Gene efficiency is held in two forms:

- `genes.eta_bar` is the raw posterior mean of the per-gene efficiency. It is
  **not** scaled by the `Inefficiency` config value.
- `genes.inefficiency` is `eta_bar * config['Inefficiency']`, the effective
  per-gene inefficiency the model actually uses. `genes.get_inefficiency()`
  returns the same thing as a gene-indexed DataFrame.

```python
obj.genes.get_inefficiency().sort_values('inefficiency').head()
```

`genes.rho_bar` is the per-gene misread density (posterior mean), also indexed by
`gene_panel`.

### Reference expression (`single_cell`)

The cell type definitions are held twice:

- `single_cell.mean_expression` is the raw reference, mean expression per gene
  (rows) per class (columns).
- `single_cell.mean_expression_adj` is the same table scaled by the
  `Inefficiency` config value, which is the version the model fits against.

### Class weights (`cellTypes`)

`cellTypes.names` are the class labels and `cellTypes.pi_bar` (also exposed as
`cellTypes.prior`) are the prior class weights, one per class, in the same order.

### Scale factors (`theta_bar`, `gamma_bar`)

The fitted scale factors, derived in [the cell scale factor](/the-model/scale-factors#theta)
and [the cell-gene scale factor](/the-model/scale-factors#gamma) sections.

- `cells.theta_bar` is `(nC, nK)`, the estimate of the cell scale factor per cell
  and class.
- `spots.gamma_bar` is `(nC, nG, nK)`, the posterior mean of gamma, the cell-gene
  scale factor, per cell, gene and class. This is a cells x genes x classes array, so it can be
  very large; index into it rather than materialising the whole thing. It is
  dropped from the pickle to keep the file small, so it is only available on
  the instance returned by `cell_type`.

## Run provenance

Every run records the version, branch and commit that produced it, in `varBayes.metadata`
and in the saved output. See [Installation](../installation.md#run-provenance).

