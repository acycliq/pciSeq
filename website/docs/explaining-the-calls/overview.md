
# Explaining the calls

pciSeq returns a probability for every cell type of every cell and for every candidate
cell of every spot. This section shows how to determine why the model arrived at those
probabilities: which genes, which prior and which neighbours determined the class of a
cell, and which terms determined the assignment of a spot to a cell or to the background.

Two methods of the fitted model answer these questions:

| question | method | page |
| --- | --- | --- |
| Why was this cell assigned this class, and not another? | [`check_cell`](../api/reference.md#check-cell) | [Why a cell is assigned its class](why-a-cell-got-its-type.md) |
| Why was this spot assigned to this cell, or to the background? | [`check_spot`](../api/reference.md#check-spot) | [Why a spot is assigned to its cell](why-a-spot-got-its-cell.md) |

The scores they show are those the model used in its last iteration, so they agree with
the probabilities in `cellData` and `geneData`.

## The fitted model

`check_cell` and `check_spot` are methods of the fitted model, a `VarBayes` object.
[`fit`](../api/reference.md#fit) does not return it, but with `save_data=True`, the
default, it is saved to `<output_path>/pciSeq/data/debug/pciSeq.pickle`:

```python
import pandas as pd

obj = pd.read_pickle('<output_path>/pciSeq/data/debug/pciSeq.pickle')
```

[`cell_type`](../api/reference.md#cell-type) returns the model as the third element of its
result. A model saved by
an earlier version of pciSeq may lack the values these methods read.

## The example data

The examples in this section use a 3D coppaFISH dataset of mouse hippocampus: 25,254
segmented cells, 3,018,812 spots of the 205 genes shared with the cell type definitions,
and 38 cell types. The data were fitted twice with identical settings except `mrf_beta`,
the strength of the spatial term: 1.5 and 0. The model from the fit with the spatial
term is `obj`, the model from the fit without it is `obj_nomrf`:

```python
import pandas as pd

obj = pd.read_pickle('espio/pciSeq/data/debug/pciSeq.pickle')              # mrf_beta = 1.5
obj_nomrf = pd.read_pickle('espio_noMRF/pciSeq/data/debug/pciSeq.pickle')  # mrf_beta = 0
```

A fitted model of this size takes several gigabytes of memory, so load one at a time if
memory is limited.
