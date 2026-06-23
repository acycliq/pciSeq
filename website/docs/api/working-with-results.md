---
title: Working with results
description: What pciSeq.fit() returns, the columns of cellData and geneData, and how to inspect the fitted model.
---

# Working with results

This page explains what comes out of a run: the two DataFrames
[`fit`](./reference#fit) returns, how their columns relate, and how to reach the
fitted model object for a deeper look.

## What `fit` returns

```python
cellData, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=ref)
```

`cellData` has one row per cell; `geneData` has one row per spot. Several columns
hold **parallel lists**: two columns of the same length whose entries line up by
position, so element `i` of one matches element `i` of the other.

## cellData

One row per cell. The background pseudo-cell is dropped, so every row is a real
cell.

| Column | Type | Meaning |
| --- | --- | --- |
| `Cell_Num` | int | Cell label. The cell's identity, matching the label image. |
| `X`, `Y`, `Z` | float | Fitted cell centroid. `Z` is present only for 3D data. |
| `ClassName` | list of str | Candidate cell classes, sorted by descending probability. |
| `Prob` | list of float | Probability of each class, lined up with `ClassName`. |
| `Genenames` | list of str | Genes assigned to the cell, sorted by descending count. |
| `CellGeneCount` | list of float | Expected count for each gene, lined up with `Genenames`. Counts are soft, so they can be fractional. |
| `spot_id` | list of list of int | The spot ids that make up each gene's count, lined up with `Genenames`. |
| `gaussian_contour` | list | The 3-sigma ellipse outline of the cell, for drawing. |
| `sphere_scale`, `sphere_rotation` | list | 3D ellipsoid drawing parameters. 3D data only. |

::: tip Coordinate units and the z-plane
`X`, `Y` and `Z` are isotropic coordinates in units of the x-voxel
(`voxel_size[0]`), not microns. pciSeq rescales every axis to the x-voxel size so
Euclidean distances are comparable in all directions. With the default isotropic
`voxel_size` of `[1, 1, 1]` the values are plain pixels and `Z` equals the
z-plane index. With anisotropic voxels, `Z` is the plane index scaled by
`voxel_size[2] / voxel_size[0]`, so to recover the (fractional) plane of a
centroid:

```python
plane = Z * voxel_size[0] / voxel_size[2]
```

Spots in `geneData` carry the integer plane directly in `plane_id`, so they need
no conversion.
:::

### How the columns relate

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
| `x`, `y`, `z` | float | Spot coordinates. `z` is present only for 3D data. |
| `plane_id` | int | The z-plane the spot sits on. |
| `neighbour` | int | The most likely parent cell label (the top of `neighbour_array`). |
| `neighbour_array` | list of int | The candidate parent cells, sorted by descending probability. The background cell `0` means a misread. |
| `neighbour_prob` | list of float | Probability of each candidate, lined up with `neighbour_array`. |
| `omp_score`, `omp_intensity` | float | Spot detection score and intensity. 3D data only. |

### How the columns relate

`neighbour_array` and `neighbour_prob` are a parallel pair: `neighbour_array[i]`
is a candidate cell and `neighbour_prob[i]` is the probability the spot belongs
to it. `neighbour` is just `neighbour_array[0]`, the winner.

The background cell, label `0`, is always one of the candidates. Its probability
is the chance the spot is a misread, so a spot with a high probability on `0` was
most likely noise.

```python
row = geneData.iloc[0]
dict(zip(row['neighbour_array'], row['neighbour_prob']))
# {458: 0.88, 12: 0.09, 0: 0.03}   # cell 458 wins; 0 is the misread chance
```

## Inspecting the fitted model

`cellData` and `geneData` are summaries. The full state of the run lives on the
fitted `VarBayes` model. See [`VarBayes`](./reference#varbayes) for the list of
attributes and methods; this section shows how to make sense of the main ones.

[`cell_type`](./reference#cell-type) returns the model directly. `fit` does not,
but with `save_data=True` (the default) it serialises the model to
`<output_path>/data/debug/pciSeq.pickle`:

```python
import pandas as pd
obj = pd.read_pickle('<output_path>/data/debug/pciSeq.pickle')
```

The model holds the run dimensions and a small object graph: `obj.cells`,
`obj.spots`, `obj.genes`, `obj.single_cell`, `obj.cellTypes`.

```python
obj.nC, obj.nS, obj.nG, obj.nK   # cells, spots, genes, classes
```

### Class probabilities (`cells.classProb`)

`cells.classProb` is a `(nC, nK)` array: one **row per cell**, one **column per
class**, and each row sums to 1. Row `0` is the background. The
columns are labelled by `cells.class_names` (the last one is `Zero`). Wrapping it
in a DataFrame makes it readable:

```python
probs = pd.DataFrame(obj.cells.classProb, columns=obj.cells.class_names)
probs.iloc[1:].idxmax(axis=1)   # most likely class for each real cell
```

`cellData`'s `ClassName` and `Prob` columns are built from this array: for each
cell the classes are sorted by descending probability and any below `0.001` are
dropped. `cells.classProb` is the raw form, keeping every class in the fixed
`class_names` order.

### Spot assignment (`spots.parent_cell_prob`)

`spots.parent_cell_prob` is `(nS, nNeighbors + 1)`: for each spot, the
probability over its nearest candidate cells. The matching cell labels are in
`spots.parent_cell_id`. The **last column is the background**, i.e. the misread
probability. `geneData`'s `neighbour_array` and `neighbour_prob` are the same
numbers sorted by descending probability; these two arrays are the raw form, with
the candidates in a fixed order.

### Genes (`genes.gene_panel`, efficiency, misread density)

`genes.gene_panel` is the master list of genes, length `nG`. Every gene-indexed
array below is in this order.

Gene efficiency comes in two forms, which are easy to confuse:

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

The single-cell reference is held twice:

- `single_cell.mean_expression` is the raw reference, mean expression per gene
  (rows) per class (columns).
- `single_cell.mean_expression_adj` is the same table scaled by the
  `Inefficiency` config value, which is the version the model fits against.

### Class weights (`cellTypes`)

`cellTypes.names` are the class labels and `cellTypes.pi_bar` (also exposed as
`cellTypes.prior`) are the prior class weights, one per class, in the same order.

### Scale factors (`theta_bar`, `gamma_bar`)

These are the fitted scale factors. They are advanced, and the maths behind them
is in [the cell scale factor](/the-model/scale-factors#theta) and
[the cell-gene scale factor](/the-model/scale-factors#gamma) pages.

- `cells.theta_bar` is `(nC, nK)`, the posterior mean cell scale factor per cell
  and class.
- `spots.gamma_bar` is `(nC, nG, nK)`, the posterior mean expression rate per
  cell, gene and class. This is a cells x genes x classes array, so it can be
  very large; index into it rather than materialising the whole thing.
