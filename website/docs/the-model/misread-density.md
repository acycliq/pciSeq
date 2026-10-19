
# Derivation: the gene-indexed misread density

A spot is assigned either to a nearby cell or, as an artefact (a misread), to the
background. A cell scores the spot by the Gaussian fit plus four expression terms; the
background scores it by the expectation (in log-space) of the misread density. Here the misread density varies
**per gene**: each gene $g$ has its own background rate $\rho_g$. The
[how it works page](../how-it-works/misread-density.md) describes the step; this page
derives the variational posterior.

## Prior

Each $\rho_g$ has a conjugate Gamma prior:

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

$r_\rho$ appears in both the shape and the rate and cancels in the prior mean, so it
changes only the concentration of the prior around $\rho_0$. The choice
$r_\rho = 1,\ \beta_\rho = 1/\rho_0$ of the source derivation is the unit-strength case.

## Likelihood

The background spots are those assigned to the "cell" $c = 0$. For gene $g$ their
log-likelihood is

$$
\log p(\mathcal{X}_{g,\,c=0} \mid \rho_g)
= \sum_{s:\, g_s = g} z_{s,0}\, \log \rho_g - \int_{\text{ROI}} \rho_g\, dx .
$$

Writing $A_{\text{total}} = \int_{\text{ROI}} dx$ for the extent of the image, which
normally coincides with the tissue section (its area in 2D, its volume in 3D), the
integral simplifies to $\rho_g A_{\text{total}}$.

### $A_{\text{total}}$ in 3D: correcting for voxel anisotropy

In 3D the ROI is a volume rather than an area, and the units it is measured in have to match
the units the rest of the score uses. They do not match by default, because confocal voxels
are not cubes: a typical `voxel_size` of `[0.28, 0.28, 0.7]` is two and a half times longer
in $z$ than in $x$ or $y$.

The position term for a spot belonging to a cell is a gaussian evaluated on
**anisotropy-scaled** coordinates, with $z$ stretched by $s_z = \texttt{voxel\_size}[2] /
\texttt{voxel\_size}[0]$ so that a micron is a micron in every direction. That gaussian is
therefore a density per unit *scaled* volume. The background term it is compared against
inside the same softmax is $\bar N_{0,g} / A_{\text{total}}$, so $A_{\text{total}}$ must be a
scaled volume too:

$$
A_{\text{total}} = W \cdot H \cdot n_{\text{planes}} \cdot s_z .
$$

Counting the ROI in raw voxels instead leaves the two sides in different units and makes the
background too strong by exactly $s_z$, which for a `[0.28, 0.28, 0.7]` voxel is a factor of
2.5. The visible symptom is spots being assigned to the background that should have gone to a
cell.

Only $z$ is corrected here. The pixels are taken to be square in $x$ and $y$, which holds for
every dataset the model has been run on. The scaling that produces the spot coordinates does
also multiply $y$ by $\texttt{voxel\_size}[1] / \texttt{voxel\_size}[0]$, so a dataset with
non-square pixels would need that factor included here as well. The default `voxel_size` is
`[1, 1, 1]`, so 2D and isotropic runs are unaffected.

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

The rate $\hat{\beta}_g = r_\rho/\rho_0 + A_{\text{total}}$ is the same for every gene: it
depends only on the prior and the extent. The gene-to-gene variation is in the shape
$\hat{r}_g = r_\rho + \bar{N}_{0,g}$, through the background count $\bar{N}_{0,g}$.

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
  \frac{\text{background spots of gene } g}{\text{extent of the image}} ,
  $$

  each gene's rate is its background spot count divided by the extent.

So $r_\rho$ interpolates between a shared constant, $\rho_0$ ($r_\rho \to \infty$), and a
per-gene empirical estimate ($r_\rho \to 0$).

In the [spot-to-cell assignment](spot-assignment.md) the background scores a spot with
$\rho_g$ alone, through $\mathbb{E}[\log\rho_g] = \psi(\hat{r}_g) - \log\hat{\beta}_g$,
with $\psi$ the digamma function. A gene with a high $\rho_g$ assigns a larger
share of each of its spots to the background.