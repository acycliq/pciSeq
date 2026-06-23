# Configuration (opts)

::: warning Auto-generated
This page is generated from the comments in `pciSeq/config.py` by
`website/gen_api.py`. Edit the comments in that file, not this page.
:::

Pass any of these as an `opts` dictionary to [`fit`](./reference#fit).
Anything you leave out keeps its default shown below.

```python
import pciSeq
opts = {'max_iter': 500, 'CellCallTolerance': 0.01}
cellData, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=ref, opts=opts)
```

### `exclude_genes`

**Default:** `[]`

list of genes to be excluded during cell-typing, e.g ['Aldoc', 'Id2'] to exclude all spots from Aldoc and Id2

### `max_iter`

**Default:** `1000`

Maximum number of loops allowed for the Variational Bayes to run

### `CellCallTolerance`

**Default:** `0.02`

Convergence achieved if assignment probabilities between two successive loops is less than the tolerance

### `rGene`

**Default:** `20`

A gamma distribution expresses the efficiency of the in-situ sequencing for each gene. It tries to capture
the ratio of the observed over the theoretical counts for a given gene. rGene controls the variance and
Inefficiency is the average of this assumed Gamma distribution

Example: If you expect 100 RNA molecules of a gene:
- With Inefficiency=0.2, you'll detect about 20 on average (20% detection rate)
- rGene=20 means this efficiency is fairly consistent between genes
- A lower rGene would mean more variance (e.g., some genes at 5% efficiency, others at 35%)

This helps account for systematic differences in detection efficiency between genes
when making cell type assignments.

### `Inefficiency`

**Default:** `0.2`

### `InsideCellBonus`

**Default:** `0`

If a spot is inside the cell boundaries this bonus will give the likelihood an extra boost
in order to make the spot more probable to get assigned to the cell than another spot positioned
outside the cell boundaries

### `mrf_beta`

**Default:** `1.0`

MRF coefficient: controls how strongly neighboring cells' class assignments
influence each other. Higher values = more spatial smoothing.

### `MisreadDensity`

**Default:** `1e-05`

MisreadDensity: Expected number of misread spots. A dictionary contains user-defined values
for gene misread densities used in the analysis.
The process to determine the misread density for each gene is as follows:

1. Compute the misread density as the number of misreads per gene divided by the area of the image.
2. Calculate the mean of these computed misread densities.
3. If the computed mean is NaN (for example, due to insufficient data), the fallback value
specified by the key 'default' in this configuration is used.
4. Finally, update the computed densities with any gene-specific overrides provided here.

The 'default' key thus sets a baseline misread density to ensure that every gene is assigned
a consistent value when the mean cannot be reliably determined.

Example configuration:
{
'default': 1e-06,  # Fallback misread density if the computed mean is NaN.
'Plp1': 0.0001,    # User-defined misread density for gene 'Plp1'.
# Additional gene-specific overrides can be added here.
}

### `cell_centroid_prior`

**Default:** `10`

A pseudo-count representing the confidence in the cell centroid estimated
by the segmentation algorithm.
- Think of `cell_centroid_prior` as if you had already seen this many imaginary
data points, all centered at the segmentation-derived centroid, before
incorporating real data.
- A small value means you have little confidence in the segmentation result,
so real data will strongly influence the estimated centroid.
- A large value means you trust the segmentation strongly, so the estimated
centroid will remain close to the initial segmentation result, even as new
data is introduced.

### `cell_cov_prior`

**Default:** `10`

A pseudo-count representing the confidence in the initial estimate of covariance,
which is modeled as a diagonal matrix where each diagonal element is the square
of the mean cell radius.
- Imagine this prior as if you had already observed this many imaginary data points
that reinforce your belief about the spread (covariance).
- A small value means you are uncertain about the initial covariance estimate,
allowing real data to significantly influence the updated covariance.
- A large value means you strongly trust the initial covariance assumption, so
real data will only cause gradual updates.

### `SpotReg`

**Default:** `0.1`

Gene detection might come with irregularities due to technical errors. A small value is introduced
here to account for these errors. It is an additive factor, applied to the single cell expression
counts when the mean counts per class and per gene are calculated.
It is like a tiny safety cushion for gene counts and adds a tiny number to all our counts to
help handle these small errors.
This is especially helpful when we see zero counts, as it prevents mathematical problems
when we're doing calculations with these numbers.

### `nNeighbors`

**Default:** `6`

By default, only the 6 nearest cells will be considered as possible parent cells for any given spot.
There is also one extra 'super-neighbor', which is always a neighbor to the spots so we can assign
the misreads to. Could be seen as the background. Hence, by default the algorithm tries examines
whether any of the 3 nearest cells is a possible parent cell to a given cell or whether the spot is
a misread

### `rSpot`

**Default:** `2`

A gamma distributed variate from Gamma(rSpot, 1) is applied to the mean expression, hence the counts
are distributed according to a Negative Binomial distribution.
The value for rSpot will control the variance/dispersion of the counts
rSpot controls how much variation we expect to see in gene counts between cells of the same type.
It's used in a Negative Binomial distribution which models gene expression.

Rules of thumb for setting rSpot:
- Default (2) is good for typical single-cell RNA data
- Lower values (0.5-1) mean high variability between cells
-> Use when you expect cells of the same type to show very different expression levels
-> Good for genes that tend to burst in expression
- Higher values (3-5) mean less variability between cells
-> Use when you expect cells of the same type to have similar expression levels
-> Good for housekeeping genes or very stable markers

Examples:
rSpot = 0.5: Counts might vary a lot (e.g., [0,5,20,100] for same cell type)
rSpot = 2.0: Moderate variation (e.g., [10,15,20,25] for same cell type)
rSpot = 5.0: More consistent counts (e.g., [17,18,19,21] for same cell type)

### `save_data`

**Default:** `True`

Boolean, if True the output will be saved as tsv files in a folder named 'pciSeq' in your system's temp dir.

### `verbose`

**Default:** `False`

Boolean. If True, turn on verbose monitoring: the per-step timings, the ELBO, and
the per-iteration diagnostic logging. None of it affects the result, and the ELBO
is expensive (several passes over the nC x nG x nK tensor), so it is off by default.

### `output_path`

**Default:** `'default'`

Set here where the results will be saved. If default then they will be saved at your system's temp folder

### `cell_radius`

**Default:** `None`

cell radius. If None then pciSeq will calc that as the mean radius across all cells.
Otherwise it will use the value provided below

### `cell_type_prior`

**Default:** `'uniform'`

cell type prior: The prior distribution on the classes. It can be 'uniform' or 'weighted'.
'uniform': weights stay fixed throughout the algorithm.
'weighted': the real class weights are updated each iteration via a Dirichlet distribution.
The Zero class weight always stays fixed (it is not part of the Dirichlet).
In both modes the initial weights come from cell_type_weights (or defaults if None).

### `cell_type_weights`

**Default:** `{'Zero': 0.5}`

cell_type_weights: A dictionary of prior probabilities for cell types.
If None: Zero = 0.5, the remaining 0.5 is split equally across real classes.
If set: specify probabilities for any subset of classes. Unspecified classes
share the remaining probability equally. Zero defaults to 0.5 if not specified.
Example: {"Zero": 0.4, "037 DG Glut": 0.1} gives Zero 40%, DG Glut 10%,
and the rest share the remaining 50% equally.

### `voxel_size`

**Default:** `[1, 1, 1]`

*******************************************************************************
Hyperparameters below added for 3D
*******************************************************************************
voxel_size: Physical size of voxels in each dimension [x, y, z].
Used to correct for anisotropic sampling in 3D data.
Example: For a microscope with:
- xy resolution of 0.147 µm/pixel
- z-step size of 0.9 µm
Use: [0.147, 0.147, 0.9]
Default: [1, 1, 1] (isotropic voxels)

### `remove_flat_cells`

**Default:** `True`

remove_flat_cells: Controls removal of cells that appear in only one z-plane.
These single-plane cells are often artifacts from segmentation, especially in 3D data.
Note: Only relevant for 3D data (multiple z-planes).
For 2D data, this setting has no effect.

### `is3D`

**Default:** `None`

Runtime attribute (automatically set during execution)

### `rTheta`

**Default:** `25.0`

### `rRho`

**Default:** `1000.0`

Shape parameter for the Gamma prior on gene-specific misread density (rho_g).
Higher values anchor rho_g closer to the MisreadDensity prior mean.
With rRho=1 the prior is weak and the data drives the estimate.

### `similarity_pairs`

**Default:** `None`

### `realtime_viewer`

**Default:** `False`

*******************************************************************************
Realtime viewer (optional visualization feature)
*******************************************************************************
Enable real-time visualization during algorithm execution

### `realtime_viewer_port`

**Default:** `5001`

Port for the realtime viewer web server

### `realtime_viewer_max_cells`

**Default:** `None`

Maximum number of cells to display (None = show all)

### `realtime_viewer_fixed_radius`

**Default:** `None`

Fixed radius for cell visualization (None = auto-scale)

### `mean_gene_counts_per_class`

**Default:** `60`

*******************************************************************************
Hyperparameters below come into action **ONLY** if single cell data are missing
*******************************************************************************

### `mean_gene_counts_per_cell`

**Default:** `30`
