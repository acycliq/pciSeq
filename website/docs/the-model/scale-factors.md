---
id: scale-factors
title: The scale factors
sidebar_label: Overview
sidebar_position: 1
---

# The scale factors

Three corrections sit between the scRNA-seq reference and the in situ experiment, each
rescaling the expected expression at a different level of detail. They are what let the
model compare a cell's observed counts against the reference on a fair footing: the
[warping the reference](../how-it-works/warping-the-reference.md) page gives the intuition,
and the pages below give the derivations.

| Factor | Indexed by | Corrects | Treated as |
| --- | --- | --- | --- |
| **[$\theta_c$](scale-theta.md)** | cell (per candidate class) | the whole cell's total count | point estimate (MAP) |
| **[$\gamma_{g,c}$](scale-gamma.md)** | gene and cell (given the class) | one gene in one cell | full Gamma (integrated out) |
| **[$\eta_g$](scale-eta.md)** | gene only (global) | one gene across the whole experiment | full Gamma |

They stack from broad to fine. $\theta_c$ is a single number for an entire cell;
$\gamma_{g,c}$ refines that down to each gene in that cell; $\eta_g$ runs the other way,
sharing one detection rate for a gene across every cell. Each absorbs the mismatch at its
own scale, and the rest is left to the others.

One of them is structurally special. The cell-gene factor $\gamma_{g,c}$ is kept as a full
Gamma random variable so it can be **integrated out**, collapsing the Poisson count model
into the Negative Binomial that the [cell-class assignment](cell-class.md) scores against.
That is why $\theta_c$ (and the latent shift $\mathbf{b}_c$ in the
[gene-gene extension](gene-gene.md)) are instead held as point estimates: keeping more than
one mixing distribution in the rate would destroy that clean collapse.

## The derivations

- **[$\theta_c$ - the cell scale factor](scale-theta.md)** - the per-cell MAP point estimate.
- **[$\gamma_{g,c}$ - the cell-gene scale factor](scale-gamma.md)** - the Poisson-Gamma
  mixture and the Negative Binomial it produces.
- **[$\eta_g$ - the in situ efficiency](scale-eta.md)** - the per-gene, reparameterised
  detection rate.
