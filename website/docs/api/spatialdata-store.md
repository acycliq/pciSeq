---
title: SpatialData store
description: The zarr store a run writes next to the tsv files, element by element, and how to query it.
---

# SpatialData store

When `save_data` is on, a run writes its results twice: as tsv and feather files, and
as a [SpatialData](https://scverse-spatialdata.readthedocs.io/) zarr store at
`<output_path>/pciSeq/data/spatialdata.zarr`. The store is for the scverse tools: napari
via napari-spatialdata, squidpy, scanpy, anndata.

This page was written against spatialdata 0.5. The elements and their columns are
produced by `pciSeq/src/core/utils/spatialdata_export.py`, which is authoritative.

## Scope

The zarr store is an optional export for tools that read the scverse formats. It
duplicates results that are already available in the tsv files and the viewer, so it
is not required to inspect a run.

Two typical uses:

**3D visualisation.** Requires
[napari-spatialdata](https://spatialdata.scverse.org/projects/napari/)
(`pip install napari-spatialdata`):

```bash
python -m napari_spatialdata view <output_path>/pciSeq/data/spatialdata.zarr
```

**Programmatic access.** See the [worked example](#worked-example) below.

## Opening the store

```python
import spatialdata as sd

sdata = sd.read_zarr('/path/to/output/pciSeq/data/spatialdata.zarr')
print(sdata)
```

```text
SpatialData object, with associated Zarr store: .../spatialdata.zarr
├── Labels
│     └── 'cell_labels': DataArray[zyx] (101, 4412, 6412)
├── Points
│     └── 'transcripts': DataFrame with shape: (<Delayed>, 10) (3D points)
└── Tables
      ├── 'cells': AnnData (26283, 207)
      └── 'reference': AnnData (38, 207)
with coordinate systems:
    ▸ 'microns', with elements:
        cell_labels (Labels), transcripts (Points)
```

`print(sdata)` lists the elements and `print(sdata.tables['cells'])` the columns and
keys of the table.

## The elements

### `cell_labels` (Labels)

The segmentation, one uint32 label image over `(z, y, x)`, or `(y, x)` for a
single-plane run. The value at a voxel is the `Cell_Num` of the cell covering it, in
the original label numbering. The element carries a transformation into the `microns`
coordinate system, so napari-spatialdata renders it at true scale next to the
transcripts.

### `transcripts` (Points)

One row per spot. Backed by parquet and lazy: call `.compute()` once to get a
pandas DataFrame.

| Column | Type | Meaning |
| --- | --- | --- |
| `x`, `y` | float32 | Spot position in pixels. |
| `z` | float32 | The z-plane index. Not the scaled `z` of `geneData`: the store keeps the plane, which is what the label image indexes. |
| `gene_name` | str | Gene name. Registered as the feature key, so `sdata.aggregate` style joins work. |
| `spot_id` | int | Spot id, same as `geneData.spot_id`. |
| `neighbour` | int | `Cell_Num` of the most likely parent cell for this spot. |
| `neighbour_prob` | float32 | Probability of that assignment. |
| `inside_cell` | int | The cell whose segmentation mask the spot falls in, as a `Cell_Num`, `0` for none. Physical containment, not an assignment. |
| `omp_score`, `omp_intensity` | float32 | Per-spot OMP diagnostics. 1.0 when the input spots had no such columns. |
| `is_hard_misread` | bool | The spot's most likely parent is the background, not a cell. |

The list columns of `geneData`, `neighbour_array` and `neighbour_prob`, are not stored,
because parquet does not store variable-length lists efficiently. The probability of the assigned cell is kept
as a scalar.

### `tables['cells']` (AnnData)

The cell typing results, one row per cell: the content of `cellData.tsv` in the anndata
layout.

| Where | What |
| --- | --- |
| `obs` index | `Cell_Num`, **as a string**. `.loc[12265]` fails, `.loc['12265']` works. |
| `obs['cell_num']` | The same number as an int. |
| `obs['region']`, `obs['class_name']`, `obs['class_prob']` | The labels element each row annotates (`'cell_labels'`), the most probable class, and its probability. |
| `obsm['spatial']` | Fitted cell centroid `(X, Y, Z)`. `X`, `Y` equal `cellData`'s; `Z` is in plane units, see the tip below. |
| `obsm['class_prob']` | The full class posterior, shape `(n_cells, n_classes)`. Columns line up with `uns['class_names']`. |
| `uns['class_names']` | Class names, the classes of the cell type definitions plus `'Zero'`, the class that expects no expression, appended last. |
| `X` | Sparse cell by gene matrix (csr) of soft counts, the same values as `cellData.CellGeneCount`. |
| `var_names` | The gene panel, matching `cellData.Genenames`. |

::: tip Z units differ from cellData
`cellData.Z` is scaled by `voxel_size[2] / voxel_size[0]` so that distances are
Euclidean. The store's `obsm['spatial']` divides it back, because the store also
holds the label image, which indexes by plane. To go from the store to the
`cellData` value: `Z * voxel_size[2] / voxel_size[0]`. With isotropic voxels the
two are the same number.
:::

The table is an annotating table (`region = 'cell_labels'`,
`instance_key = 'cell_num'`), so scanpy and napari-spatialdata can join it
against the segmentation automatically.

### `tables['reference']` (AnnData)

The cell type definitions the run was scored against, transposed to classes as rows and
genes as columns: `X[i, j]` is the mean expression of gene `j` in class `i`.

## Provenance

`sdata.attrs['pciseq']` records what produced the store: the code (`branch`,
`commit`, `commit_date`, `version`), the environment (`python_version`, `os`,
`package_versions`), when the run was created (`created_at`), the
`spatialdata_version`, and the resolved run `config`, without `label_map`.

## Mapping from the tsv columns

| tsv column | in the store |
| --- | --- |
| `cellData.Cell_Num` | `tables['cells'].obs` index, `obs['cell_num']`, label image values |
| `cellData.X`, `Y`, `Z` | `obsm['spatial']` (Z in plane units) |
| `cellData.ClassName`, `Prob` | `obsm['class_prob']` + `uns['class_names']`; argmax also flattened into `obs` |
| `cellData.Genenames`, `CellGeneCount` | `X` + `var_names` |
| `geneData` columns | `transcripts` points, see the table above |

## Worked example

One cell, looked up by `Cell_Num`, across the whole store:

```python
import numpy as np

adata = sdata.tables['cells']
cid = 12265
pos = int(np.flatnonzero(adata.obs.index.values == str(cid))[0])

obs = adata.obs.iloc[pos]        # class_name, class_prob, ...
xyz = adata.obsm['spatial'][pos] # fitted cell centroid, Z in plane units
prob = adata.obsm['class_prob'][pos]
names = np.asarray(adata.uns['class_names'])
top3 = np.argsort(-prob)[:3]
[(names[i], round(float(prob[i]), 3)) for i in top3]
# [('037 DG Glut', 1.0), ('038 DG-PIR Ex IMN', 0.0), ('017 CA3 Glut', 0.0)]

counts = np.asarray(adata.X[pos].todense()).ravel()
adata.var_names[np.argmax(counts)], counts.max()   # top gene in this cell
# ('Sema5a', 6.0)
```

Spots whose most likely parent is this cell, and the label image check at the
fitted cell centroid:

```python
tr = sdata['transcripts'].compute()
tr[tr.neighbour == cid]

x, y, z = adata.obsm['spatial'][pos]
int(sdata['cell_labels'][int(z), int(y), int(x)])   # usually cid; the fitted centroid need not lie inside the mask
# single-plane run: sdata['cell_labels'][int(y), int(x)]
```

::: warning Two ways to count a cell's spots
`len(tr[tr.neighbour == cid])` (hard argmax assignment) and
`adata.X[pos].sum()` (the model's soft counts) are close but not identical. They
answer different questions. The tsv files and the store agree on both.
:::
