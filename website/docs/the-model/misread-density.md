---
id: misread-density
title: Derivation - gene-indexed misread density
sidebar_label: Misread density
sidebar_position: 2
---

# Derivation: the gene-indexed misread density

A spot is either claimed by a nearby cell or treated as an artefact (a misread) and
assigned to the background. A cell claims a spot through its spatial term combined with
how well the spot fits the cell's likely class; the background claims it through the mean
of the misread density. Here the misread density varies **per gene**: each gene $g$ has
its own background rate $\rho_g$. The [how it works page](../how-it-works/misread-density.md)
gives the intuition; this page derives the variational posterior.

## Prior

We place a conjugate Gamma prior on each $\rho_g$:

$$
p(\rho_g) = \mathrm{Gamma}(\rho_g;\, r_\rho, \beta_\rho)
= \frac{\beta_\rho^{\,r_\rho}}{\Gamma(r_\rho)}\,
  \rho_g^{\,r_\rho - 1}\, e^{-\beta_\rho \rho_g},
$$

where $r_\rho$ is the shape and $\beta_\rho$ the rate, giving prior mean
$\mathbb{E}[\rho_g] = r_\rho/\beta_\rho$.

Two of these symbols map directly onto configuration settings in the code:

- $r_\rho$ is the **`rRho`** setting (the prior strength);
- $\rho_0$ is the **`MisreadDensity`** setting (the prior mean misread density).

The implementation parameterises the prior so that its **mean is fixed** at $\rho_0$, while
$r_\rho$ controls how strongly that mean is held. The shape is `rRho` directly, and the
rate is chosen to hold the mean at $\rho_0$:

$$
\beta_\rho = \frac{r_\rho}{\rho_0},
\qquad\text{so that}\qquad
\mathbb{E}[\rho_g] = \frac{r_\rho}{\beta_\rho} = \rho_0 \quad\text{for any } r_\rho .
$$

Because $r_\rho$ appears in both the shape and the rate, it cancels in the prior mean:
changing $r_\rho$ moves no probability mass off $\rho_0$, it only changes the
**concentration** of the prior around it. This is what makes $r_\rho$ a clean "prior
strength" dial, as the posterior below makes explicit. (The simpler choice
$r_\rho = 1,\ \beta_\rho = 1/\rho_0$ in the source derivation is the special case of this at
unit strength.)

## Likelihood

The background spots are those assigned to the "cell" $c = 0$. For gene $g$ their
log-likelihood is

$$
\log p(\mathcal{X}_{g,\,c=0} \mid \rho_g)
= \sum_{s:\, g_s = g} z_{s,0}\, \log \rho_g - \int_{\text{ROI}} \rho_g\, dx .
$$

Writing $A_{\text{total}} = \int_{\text{ROI}} dx$ for the total area of the tissue
section, the integral simplifies to $\rho_g A_{\text{total}}$.

## Variational update

Under the mean-field approximation, the update for $q(\rho_g)$ is the expected log-joint
over all other latent variables, keeping only the terms that depend on $\rho_g$:

$$
\log q(\rho_g)
= \mathbb{E}_{z,\zeta,\gamma,\eta,\theta}\!\left[
  \sum_{s:\, g_s=g} z_{s,0}\, \log \rho_g - \rho_g A_{\text{total}} + \log p(\rho_g)
\right] + \text{const}.
$$

Taking the expectation replaces the spot indicators by their expected counts. Let

$$
\bar{N}_{0,g} = \sum_{s:\, g_s = g} q(z_{s,0} = 1)
$$

be the expected number of spots of gene $g$ assigned to the background. Substituting the
Gamma prior $\log p(\rho_g) = (r_\rho - 1)\log\rho_g - \beta_\rho \rho_g$:

$$
\log q(\rho_g)
= \bar{N}_{0,g}\, \log \rho_g - \rho_g A_{\text{total}}
  + (r_\rho - 1)\log \rho_g - \beta_\rho \rho_g + \text{const}.
$$

Grouping the $\log\rho_g$ and $\rho_g$ terms:

$$
\log q(\rho_g)
= (\bar{N}_{0,g} + r_\rho - 1)\, \log \rho_g
  - (A_{\text{total}} + \beta_\rho)\, \rho_g + \text{const}.
$$

This is the log-density of a Gamma distribution. Therefore

$$
\boxed{\;
q(\rho_g) = \mathrm{Gamma}(\rho_g;\, \hat{r}_g, \hat{\beta}_g)
\;}
$$

with updated parameters. Substituting the implementation's prior rate
$\beta_\rho = r_\rho/\rho_0$:

$$
\hat{r}_g = r_\rho + \bar{N}_{0,g},
\qquad
\hat{\beta}_g = \frac{r_\rho}{\rho_0} + A_{\text{total}} .
$$

## Reading the result

The posterior mean is

$$
\mathbb{E}[\rho_g]
= \frac{\hat{r}_g}{\hat{\beta}_g}
= \frac{r_\rho + \bar{N}_{0,g}}{\dfrac{r_\rho}{\rho_0} + A_{\text{total}}} .
$$

Notice that the rate $\hat{\beta}_g = r_\rho/\rho_0 + A_{\text{total}}$ is the **same for
every gene** - it depends only on the prior and the tissue area, not on $g$. All the
gene-to-gene variation lives in the shape $\hat{r}_g = r_\rho + \bar{N}_{0,g}$, through the
background count $\bar{N}_{0,g}$.

### The role of $r_\rho$ (the prior strength)

The parameter $r_\rho$ (the `rRho` setting) acts as a pseudo-count that sets how readily the
data move $\rho_g$ away from the prior mean $\rho_0$ (the `MisreadDensity` setting):

- **Large $r_\rho$ - the prior dominates and the misread density stops being gene-specific.**
  The shape $r_\rho + \bar{N}_{0,g} \approx r_\rho$ for every gene, because the background
  counts are swamped by the large pseudo-count, so every gene's posterior mean sits at
  $\rho_0$ - the single global constant of the original model.

- **Small $r_\rho$ - the misread density is gene-specific and fully data driven.** With the
  prior contributing very weakly, the shape
  $r_\rho + \bar{N}_{0,g} \approx \bar{N}_{0,g}$ and the rate
  $r_\rho/\rho_0 + A_{\text{total}} \approx A_{\text{total}}$ (e.g. at $r_\rho = 1$), so

  $$
  \mathbb{E}[\rho_g] \approx \frac{\bar{N}_{0,g}}{A_{\text{total}}}
  \;=\;
  \frac{\text{background spots of gene } g}{\text{tissue area}} ,
  $$

  each gene's noise floor is just its background spot count divided by the tissue area.

So $r_\rho$ interpolates between a shared constant ($r_\rho \to \infty$) and a per-gene
empirical estimate ($r_\rho \to 0$), with $\rho_0$ as the anchor in both limits.

The misread density is the background's means of claiming a spot. Every spot is contested
in the [spot-to-cell assignment](../how-it-works/spots-to-cells.md): it is weighed against
each neighbouring cell and against the background, and assigned to whichever makes the
strongest case. A cell makes its case with several quantitative and qualitative factors -
how close the spot is, and how well its gene fits the cell's likely type. The background
has only one: the misread density $\rho_g$. A noisy gene has a high $\rho_g$, so the
background presses a stronger claim, and a genuine spot of that gene must beat that higher
bar to be won by a cell. This is the value that enters the assignment as the background
score, through $\mathbb{E}[\log\rho_g] = \psi(\hat{r}_g) - \log\hat{\beta}_g$ (with $\psi$
the digamma function).