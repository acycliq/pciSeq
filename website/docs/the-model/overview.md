---
id: overview
title: The model - formal definition
sidebar_label: Overview
sidebar_position: 1
---

# The model: formal definition

The [how it works](../how-it-works/overview.md) section explains the algorithm one block
at a time, in words. This section states the same model formally and derives every update
equation. It follows *The extended pciSeq model* (v0.3), which is the version the
`dev_3d` code implements. The model builds on the original construction of Qian et al.
(2020); the corrections we apply to that paper are listed in the [errata](errata.md).

## Notation

| Symbol | Meaning |
| --- | --- |
| $x_s,\ g_s$ | location and decoded gene label of RNA spot $s$ |
| $z_{s,c}$ | indicator assigning spot $s$ to cell $c$ |
| $\zeta_{c,k}$ | indicator assigning cell $c$ to class $k$ |
| $\mu_{g,k}$ | mean expression of gene $g$ in class $k$ (from scRNA-seq) |
| $\gamma_{g,c}$ | per-gene, per-cell scale factor |
| $\eta_g$ | in situ detection efficiency of gene $g$ |
| $\eta_0$ | global prior mean of detection efficiency (typically 0.2) |
| $\theta_c$ | per-cell scale factor (extension) |
| $\rho_g$ | per-gene background (misread) density (extension) |
| $A_c=\int e^{-D_c(x)}\,dx$ | normalised area of cell $c$ |
| $D_c(x)$ | distance from point $x$ to cell $c$ |
| $\pi_k$ | prior probability of class $k$ |
| $\mathcal{N}_c$ | the set of nearest neighbours of cell $c$ |
| $\beta$ | strength of the spatial (MRF) coupling |

## The generative model

RNA spots of gene $g$ arising from cell $c$ are modelled as a spatial Poisson process
with intensity

$$
\lambda_{g,c}(x) = \mu_{g,k(c)}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g .
$$

The log-likelihood of a Poisson process has the general form
$-\int \lambda(x)\,dx + \sum_s \log \lambda(x_s)$: the integral charges a cost
proportional to the expected total count, and each observed spot contributes its
log-intensity. Summing over genes, cells, and classes and adding the priors on
$\gamma$, $\eta$, and the class assignments $\zeta$, the original log-joint of
Qian et al. (2020) is

$$
\begin{aligned}
\log p(x, g, z, \zeta, \gamma, \eta) =
& - \sum_{g,c,k} \zeta_{c,k}\, \mu_{g,k}\, A_c\, \gamma_{g,c}\, \eta_g \\
& + \sum_{s,c,k} z_{s,c}\, \zeta_{c,k}\,
   \log\!\big[\, \mu_{g_s,k}\, e^{-D_c(x_s)}\, \gamma_{g_s,c}\, \eta_{g_s} \big] \\
& + \sum_{g,c} \log p(\gamma_{g,c}) + \sum_g \log p(\eta_g)
   + \sum_{c,k} \zeta_{c,k}\, \log \pi_k .
\end{aligned}
$$

The extended model adds two ingredients: a **per-cell scale factor** $\theta_c$ inside the
intensity, and a **Markov Random Field** prior that rewards neighbouring cells for sharing
a class. The log-joint becomes

$$
\begin{aligned}
\log p(x, g, z, \zeta, \gamma, \eta, \theta) =
& - \sum_{g,c,k} \zeta_{c,k}\, \theta_c\, \mu_{g,k}\, A_c\, \gamma_{g,c}\, \eta_g \\
& + \sum_{s,c,k} z_{s,c}\, \zeta_{c,k}\,
   \log\!\big[\, \theta_c\, \mu_{g_s,k}\, e^{-D_c(x_s)}\, \gamma_{g_s,c}\, \eta_{g_s} \big] \\
& + \sum_{g,c} \log p(\gamma_{g,c}) + \sum_g \log p(\eta_g) + \sum_c \log p(\theta_c) \\
& + \sum_{c,k} \zeta_{c,k}\, \log \pi_k
   + \beta \sum_{c,k} \sum_{c' \in \mathcal{N}_c} \mathbf{1}(\zeta_{c,k} = \zeta_{c',k} = 1) .
\end{aligned}
$$

The final term is the MRF: since $\zeta_{c,k}\in\{0,1\}$, the indicator
$\mathbf{1}(\zeta_{c,k}=\zeta_{c',k}=1)$ is simply the product $\zeta_{c,k}\zeta_{c',k}$,
and we use that product form in all derivations. The background (misread) density is also
promoted from a single global constant to a per-gene quantity $\rho_g$, derived
[on its own page](misread-density.md).

## The variational approximation

The posterior is intractable, so we approximate it by **variational inference**: we pick
the member of a tractable, factorised family that is closest to the true posterior in
Kullback-Leibler divergence, and fit it by **coordinate ascent** (CAVI), updating one
factor at a time. Because $\gamma$ is meant to depend on $\theta$ and on the cell's class,
those three are bundled into one structured factor:

$$
p(z, \zeta, \gamma, \eta, \theta \mid x, g)
\approx q(\gamma \mid \zeta, \theta)\, q(\theta \mid \zeta)\, q(\zeta)\, q(z)\, q(\eta) .
$$

One factor needs special handling. The per-cell scale $\theta_c$ enters the intensity
multiplicatively with $\gamma_{g,c}$; treating both as full random variables would make
the marginalisation intractable and destroy the Negative Binomial likelihood that drives
cell typing. We therefore restrict $q(\theta_c)$ to a **point estimate** (a Dirac delta),
which keeps $\theta_c$ constant during the update for $\gamma_{g,c}$. Mixing a Dirac factor
with full variational factors makes the scheme **Variational EM** rather than pure
variational Bayes; the [self-consistency appendix](appendix-self-consistency.md) records
why this is sound.

## What these pages derive

- **[Misread density $\rho_g$](misread-density.md)** - the per-gene background rate.
- **[Cell scale factor $\theta_c$](cell-scale-theta.md)** - the per-cell point estimate.
- **[Scale factors $\gamma_{g,c}$ and $\eta_g$](scale-factors.md)** - the per-gene-per-cell
  and per-gene corrections.
- **[Cell-class assignment $q(\zeta)$](cell-class.md)** - the Negative Binomial likelihood,
  the class prior, and the MRF spatial term.
- **[Spot-to-cell assignment $q(z)$](spot-assignment.md)** - which cell (or the background)
  each spot is assigned to.
- **[Gene-gene dependence](gene-gene.md)** - the latent shift $\mathbf{b}_c$ and its
  class covariance $\Sigma_k$.

The [errata](errata.md) lists the corrections to Qian et al. (2020), and the
[self-consistency appendix](appendix-self-consistency.md) records why the Dirac
construction is sound.

## Not yet documented

Two parts of the model are developed but not yet written up here, matching the stubs in
the source document:

- a **Dirichlet prior** on the class probabilities $\pi$;
- **cell positions and shapes as a Gaussian mixture**. The current implementation assumes
  fixed, spherical shapes; the Gaussian-mixture extension exists but is disabled in the
  Python code.