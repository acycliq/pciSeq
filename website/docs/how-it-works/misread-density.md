
# 1. Estimating the misread density

A misread is a spot that is attributed to the background rather than to a cell. Such
spots arise from imaging and decoding artefacts, and some survive the quality filters
of the spot calling.

The misread density is the rate, per gene and per unit of imaged area or volume, at
which misreads occur. Its expected log is the background's score, against which each
candidate cell's score is compared.

## The background process

The background is a spatially uniform process: noise spots of a gene occur at a
constant rate everywhere in the section. In [spot assignment](spots-to-cells.md) this
rate gives the score of the background option, the same at every location. The higher a
cell's score relative to the background's, the larger the share of the spot's probability
it receives.

## Per gene, not one global rate

The original pciSeq used one background rate for all genes. This version estimates a
rate per gene, since genes differ in how noisy they are: for a noisy gene the background
option receives a larger share of each spot's probability than for a clean one.

The rate is updated on every iteration from the spots currently attributed to the
background:

$$
\text{background rate of gene } g \;\approx\;
\frac{\text{background spots of gene } g}{\text{extent of the ROI}}
$$

With a conjugate Gamma prior the rate has a Gamma posterior whose mean is this ratio,
regularised by the prior:

$$
q(\rho_g) = \mathrm{Gamma}\Big(\rho_g;\; r_\rho + \bar{N}_{0,g},\;\; \tfrac{r_\rho}{\rho_0} + A_{\text{total}}\Big),
$$

where $\bar{N}_{0,g}$ is the expected number of background spots of gene $g$,
$A_{\text{total}}$ is the extent of the region of interest, $\rho_0$ is the prior mean
misread density, and $r_\rho$ sets how strongly that prior is held.

In 2D the extent is an area. In 3D it is a volume, measured in the same units as the
cell shapes it is compared against, so the voxel anisotropy has to be corrected for. An
extent that is too small makes the background too strong and causes spots that belong
to cells to be assigned to the background. The [derivation](../the-model/misread-density.md) covers the
correction and the role of $r_\rho$.
