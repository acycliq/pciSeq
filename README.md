# pciSeq: Probabilistic Cell typing by In situ Sequencing

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://acycliq.github.io/pciSeq_3d/)

pciSeq takes the raw output of an in situ sequencing experiment and works out two
things at the same time: which cell each RNA spot belongs to, and what cell type each
cell is. It does this with a probabilistic model that leans on a separate single-cell
RNA-seq reference, so both answers come out as probabilities rather than hard labels.

## Documentation

Full docs, including a plain-language walkthrough of how the algorithm works, are at
**[docs](https://acycliq.github.io/pciSeq_3d/)**.

[//]: # (- [What is pciSeq?]&#40;https://acycliq.github.io/pciSeq_3d/&#41;)

[//]: # (- [How it works &#40;the variational loop&#41;]&#40;https://acycliq.github.io/pciSeq_3d/how-it-works/overview&#41;)

## Install

```bash
pip install git+https://github.com/acycliq/pciSeq_3d.git@dev_3d
```

## Quick start

```python
import pciSeq

# spots: DataFrame with columns gene_name, x, y (and z_plane for 3D)
# coo: segmentation label image(s) as scipy.sparse coo_matrix
# scRNAseq: single-cell reference, genes as index, cell types as columns
cellData, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=scRNAseq)
```

## Related

- [pciSeq Viewer](https://github.com/acycliq/pciSeq_viewer) - desktop and web app for
  exploring pciSeq results.