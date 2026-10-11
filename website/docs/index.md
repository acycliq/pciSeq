# pciSeq

pciSeq (probabilistic cell typing by in situ sequencing) assigns each RNA spot of an
imaging-based spatial transcriptomics experiment to a cell, and each cell to a cell type.
The two assignments are estimated jointly: the type of a cell depends on the spots inside
it, and the cell of a spot depends on the types of the cells around it. Both are returned
as probabilities.

This version extends the method of
[Qian et al. (2020)](https://doi.org/10.1038/s41592-019-0631-4) in several ways:

- **3D data.** Spots and segmentation may be volumes, with the voxel size taken into
  account.
- **A scale factor per cell.** It absorbs differences in total count between cells, from
  cell size or capture efficiency, so two cells of the same class are recognised as such
  even when one holds far fewer spots than the other. See
  [scale factors](the-model/scale-factors.md).
- **A background density per gene, estimated from the data.** In the original it is one
  fixed value. See [misread density](the-model/misread-density.md).
- **A spatial term.** The classes of the neighbouring cells can be taken into account
  when deciding the class of a cell. See
  [cell-class assignment](the-model/cell-class.md).
- **A class prior estimated from the data**, optional. Classes that are common in the
  tissue get a larger prior. See [class prior](the-model/class-prior.md).
- **Cell geometry.** The probability that a spot comes from a cell falls with its distance
  from the cell, measured by a Gaussian. In the original the Gaussian is the same for every
  cell: centred on the segmentation centroid, with the mean cell radius. Here its centre
  and its spread are fitted per cell from the spots: it can become an ellipsoid, elongated
  along the spots, and its centre can shift where the segmentation did not capture the cell
  well. See [cell geometry](the-model/gaussian.md).

## Inputs

- **Spots.** A table of detected RNA spots with gene identity and position (`x`, `y`, and
  `z_plane` for 3D data).
- **Segmentation.** A label image giving the cell each pixel belongs to, typically from a
  DAPI nuclear stain.
- **Cell type definitions.** Mean expression per gene for each cell type, for example
  from a scRNA-seq experiment. These are the reference profiles the cells are scored against.

## Outputs

- For every cell, a probability distribution over the cell types.
- For every spot, a probability distribution over its candidate parent cells and the
  background.

`fit` returns both as pandas DataFrames. If `save_data` is `True`, they are also written
to `output_path` as tsv and feather files and as a SpatialData zarr store, and the fitted
model is pickled. See [Working with results](api/working-with-results.md).

## Probabilistic output

Assignments are probabilities rather than labels, because segmentation boundaries are
imprecise, detection is imperfect and a fraction of the spots are noise. A spot on the
boundary of two cells receives probability on both; a spot that no cell explains is
assigned to the background. Cell type assignments are probabilities in the same way.

A cell's count of a gene is the sum of the assignment probabilities of that gene's spots,
so it is an expected count rather than an integer. Four spots of Plp1 assigned to the cell
with probability 1 and one with probability 0.3 give a count of 4.3.

## Usage

[Install](installation.md), then see [Running pciSeq](running-pciseq.md). The
algorithm is described in [How it works](how-it-works/overview.md).
[Explaining the calls](explaining-the-calls/overview.md) shows how to find out why the
model assigned a cell to its type and a spot to its cell.
