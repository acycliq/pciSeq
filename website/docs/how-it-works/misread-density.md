---
id: misread-density
title: Estimating the misread density
sidebar_label: 1. Misread density
sidebar_position: 2
---

# Block 1: Estimating the misread density

Not every detected spot corresponds to a genuine transcript. Some reads arise from
technical artefacts, optical crosstalk, or decoding errors. The misread density
quantifies this background, so that each spot can be weighed against the possibility
that it is noise rather than signal.

## The idea in one line

For each gene, estimate the rate at which its spots are produced by background noise
spread uniformly across the tissue, and use that rate as the baseline a genuine spot
must exceed.

## How it works

The background is modelled as a spatially uniform process: noise spots of a given gene
occur at a constant rate everywhere in the section, independent of location. A gene with
a high background rate therefore produces spurious spots throughout the tissue, including
regions far from any cell, whereas a gene with a low rate does so rarely.

When a spot is assigned later, in [block 4](spots-to-cells.md), this rate sets the score
of the background option. The background competes with the candidate cells as a constant,
location-independent alternative: a spot is assigned to a cell only if that cell explains
it better than the background would.

## Per gene, not one global rate

The original pciSeq used a single background rate for all genes. This version estimates a
separate rate for each gene, since genes differ in how noisy they are. With a per-gene
rate, spots of a noisy gene must exceed a higher background level, while spots of a clean
gene need only exceed a lower one.

The rate is updated on every iteration from the spots currently attributed to the
background:

$$
\text{background rate of gene } g \;\approx\;
\frac{\text{background spots of gene } g}{\text{tissue area}}
$$

This ratio of an observed count to the extent over which it is spread is the first
instance of a form that recurs throughout pciSeq: an estimate expressed as an observed
quantity divided by an expected one. The same structure underlies the scaling factors in
the next block.

Formally, with a conjugate Gamma prior the background rate has a Gamma posterior whose
mean is exactly this ratio (regularised by the prior):

$$
q(\rho_g) = \mathrm{Gamma}\Big(\rho_g;\; r_\rho + \bar{N}_{0,g},\;\; \tfrac{r_\rho}{\rho_0} + A_{\text{total}}\Big),
$$

where $\bar{N}_{0,g}$ is the expected number of background spots of gene $g$,
$A_{\text{total}}$ is the tissue area, $\rho_0$ is the prior mean misread density, and
$r_\rho$ sets how strongly that prior is held. The
[full derivation](../the-model/misread-density.md), including the role of $r_\rho$, is in
the model section.

## What feeds in and what comes out

- **Feeds in:** the spots attributed to the background on the previous iteration.
- **Comes out:** a per-gene background rate, used by [block 4](spots-to-cells.md) as the
  background option against which each spot is compared.