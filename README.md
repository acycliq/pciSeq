<div align="center">

# pciSeq

**Probabilistic cell typing for imaging-based spatial transcriptomics, in 3D**

[![Documentation](https://img.shields.io/badge/docs-acycliq.github.io-0a7d5a)](https://acycliq.github.io/pciSeq_3d/)
[![Tests](https://github.com/acycliq/pciSeq_3d/actions/workflows/pytest.yaml/badge.svg?branch=dev_3d)](https://github.com/acycliq/pciSeq_3d/actions/workflows/pytest.yaml)
[![Python](https://img.shields.io/badge/python-3.10%2B-3776ab)](#installation)
[![Licence](https://img.shields.io/badge/licence-MIT-lightgrey)](LICENSE)

[Documentation](https://acycliq.github.io/pciSeq_3d/) ·
[Installation](#installation) ·
[Quick start](#quick-start)

</div>

pciSeq assigns each RNA spot of an imaging-based spatial transcriptomics experiment to a
cell, and each cell to a cell type. The two assignments are estimated jointly, in 3D,
and both come back as probabilities.

<p align="center">
<img src="https://github.com/user-attachments/assets/dfcd8438-e637-4e34-ae73-4e8c445a982b" width="100%" alt="The live viewer during a fit">
</p>

<p align="center"><sub><em>The live viewer during a fit, one frame per iteration. Cells are
coloured by the class they are currently assigned to, and the chart at the bottom right
tracks convergence.</em></sub></p>

## Installation

```bash
pip install git+https://github.com/acycliq/pciSeq_3d.git@dev_3d
```

Python 3.10 or newer. See
[installation](https://acycliq.github.io/pciSeq_3d/installation) for the optional `mcp`
extra.

## Quick start

```python
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import pciSeq

spots = pd.read_csv("spots.csv")                             # gene_name, x, y, z_plane
coo = [coo_matrix(plane) for plane in np.load("masks.npy")]  # one label matrix per plane
scRNAseq = pd.read_csv("scRNAseq.csv", index_col=0)          # genes x cell types

opts = {
    "voxel_size": [0.28, 0.28, 0.7],
    "output_path": "out/run1",
    "realtime_viewer": True,     # watch it run at http://127.0.0.1:5001
}

cellData, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=scRNAseq, opts=opts)
```

`cellData` holds, for every cell, a probability distribution over the cell types and its
expected gene counts. `geneData` holds, for every spot, a probability distribution over
its candidate parent cells and the background. Both are pandas DataFrames, and both are
also written to the output folder, as tsv and as a
[SpatialData](https://spatialdata.scverse.org) zarr store.

pciSeq saves with its results the terms from which every assignment was computed. An AI
agent, or the methods `check_cell` and `check_spot`, read them to explain why a cell was
assigned to its class or why a spot was assigned to a cell.

## Documentation

| Section | Description |
| --- | --- |
| [Running pciSeq](https://acycliq.github.io/pciSeq_3d/running-pciseq) | from Python, from the command line, the settings and how to tune them |
| [How it works](https://acycliq.github.io/pciSeq_3d/how-it-works/overview) | the building blocks of the algorithm |
| [Explaining the calls](https://acycliq.github.io/pciSeq_3d/explaining-the-calls/overview) | `check_cell` and `check_spot`, followed through on one cell and one spot |
| [Accessing the AI agent](https://acycliq.github.io/pciSeq_3d/accessing-the-ai-agent) | the AI agent: answers questions about the results or about how pciSeq works |
| [Viewer](https://acycliq.github.io/pciSeq_3d/viewer/overview) | pciSeq Viewer, the desktop application to visualise the results |
| [The model](https://acycliq.github.io/pciSeq_3d/the-model/overview) | the Bayesian model, stated formally, and the derivation of its update equations |
| [API](https://acycliq.github.io/pciSeq_3d/api/reference) | the reference, and a [code map](https://acycliq.github.io/pciSeq_3d/api/code-map) from each model quantity to the line that computes it |

## Licence

MIT, see [LICENSE](LICENSE).
