
# The model: formal definition

The [how it works](../how-it-works/overview.md) section explains the algorithm one piece
at a time, in words. This section states the same model formally and derives every update
equation. It follows *The extended pciSeq model* (v0.3), which is the version the
`dev_3d` code implements. The model builds on the original construction of Qian et al.
(2020).

## Notation

| Symbol | Meaning |
| --- | --- |
| $x_s,\ g_s$ | location and decoded gene label of RNA spot $s$ |
| $z_{s,c}$ | indicator assigning spot $s$ to cell $c$ |
| $\zeta_{c,k}$ | indicator assigning cell $c$ to class $k$ |
| $\mu_{g,k}$ | mean expression of gene $g$ in class $k$ (from the cell type definitions) |
| $\gamma_{g,c}$ | per-gene, per-cell scale factor |
| $\eta_g$ | in situ detection efficiency of gene $g$ |
| $\eta_0$ | global prior mean of detection efficiency (typically 0.2) |
| $\theta_c$ | per-cell scale factor (extension) |
| $\rho_g$ | per-gene background (misread) density (extension) |
| $A_c=\int e^{-D_c(x)}\,dx$ | a per-cell normalising factor tied to the inside-cell bonus. It equals $1$ when the bonus is off |
| $D_c(x)$ | negative log density of point $x$ under the Gaussian of cell $c$ |
| $m_c,\ \Sigma_c$ | centre and covariance of that Gaussian (extension) |
| $d$ | number of spatial dimensions, $2$ or $3$ |
| $\kappa_0,\ \nu_0$ | pseudo-count of the prior on $m_c$ and degrees of freedom of the prior on $\Sigma_c$ (extension) |
| $\pi_k$ | prior probability of class $k$; fixed, or a random variable under `weighted` (extension) |
| $\delta_k$ | Dirichlet pseudo-count of the prior on $\pi_k$ (extension) |
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

The extended model adds four ingredients: a **per-cell scale factor** $\theta_c$ inside the
intensity, a **Markov Random Field** prior that rewards neighbouring cells for sharing
a class, a prior on each **cell's geometry**, covered
[below](#cell-geometry), and a prior on the **class probabilities** $\pi$, covered
[below](#class-prior). With the first two the log-joint becomes

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
\tag{1}
$$

The final term is the MRF: since $\zeta_{c,k}\in\{0,1\}$, the indicator
$\mathbf{1}(\zeta_{c,k}=\zeta_{c',k}=1)$ is simply the product $\zeta_{c,k}\zeta_{c',k}$,
and the product form is used in all derivations.

### The background

So far every spot is taken to come from a cell: the sums over $c$ run over the real cells,
$c = 1, \dots, C$. Some spots come from no cell. They are misreads, and the model gives them
a component of their own, the background, $c = 0$.

The background is not a cell. It has no class, no position and no scale factors, so $\zeta$,
$\mu$, $D_c$, $\gamma$, $\theta$ and $\eta$ do not apply to it. Misreads of gene $g$ fall
over the whole region at a constant density, $\rho_g$. In the original model that density is
one fixed number. Here it is a random variable, one per gene, with a prior $p(\rho_g)$, and
it is estimated along with everything else.

With the background included, the intensity has two branches:

$$
\lambda_{g,c}(x) =
\begin{cases}
\theta_c\, \mu_{g,k(c)}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g , & c = 1, \dots, C \quad \text{(a cell)} \\[4pt]
\rho_g , & c = 0 \quad \text{(the background)}
\end{cases}
$$

Both are Poisson processes, so the general form of the log-likelihood,
$-\int \lambda(x)\,dx + \sum_s \log \lambda(x_s)$, applies to both. The log-joint of the top
branch is given by equation (1). For the bottom branch, write $\mathcal{X}_{g,0}$ for the
spots of gene $g$ assigned to the background. Conditional on $\rho_g$, their log-likelihood
is

$$
\log p(\mathcal{X}_{g,0} \mid \rho_g)
= - \int_{\text{ROI}} \rho_g\, dx + \sum_{s:\, g_s = g} z_{s,0}\, \log \rho_g
= - \rho_g\, A_{\text{total}} + \sum_{s:\, g_s = g} z_{s,0}\, \log \rho_g ,
$$

where $A_{\text{total}} = \int_{\text{ROI}} dx$ is the volume of the region and
$z_{s,0} = 1$ marks a spot assigned to the background. The first part is the cost of the
misreads the density predicts over the region; the second is the log-intensity of each
spot put down to the background.

Adding the prior $\log p(\rho_g)$ and summing over genes turns this into a joint. With
$\mathcal{X}_0$ for all the background spots, the log-joint of the bottom branch is

$$
\log p(\mathcal{X}_0, \rho)
= \sum_g \Big[ - \rho_g\, A_{\text{total}} + \sum_{s:\, g_s = g} z_{s,0}\, \log \rho_g
+ \log p(\rho_g) \Big] .
\tag{2}
$$

The posterior of $\rho_g$ is derived [on its own page](misread-density.md).

### Cell geometry

The position factor $e^{-D_c(x)}$ is a Gaussian density, $\mathcal{N}(x \mid m_c, \Sigma_c)$,
so $A_c = 1$ when the bonus is off. In the original model its centre $m_c$ is the
segmentation centroid and its covariance $\Sigma_c$ is $r^2 I$, with $r$ the mean cell
radius, both fixed. The extended model makes them random variables with a conjugate
Normal-Inverse-Wishart prior,

$$
p(m_c, \Sigma_c)
= \mathcal{N}\!\big(m_c \mid m_{0,c},\, \Sigma_c / \kappa_0\big)\;
  \mathcal{IW}\!\big(\Sigma_c \mid (\nu_0 - d - 1)\, r^2 I,\, \nu_0\big),
$$

centred on the segmentation centroid $m_{0,c}$ and on $r^2 I$. The centroid is worth
$\kappa_0$ spots and the shape $\nu_0 - d - 1$. This adds $\sum_c \log p(m_c, \Sigma_c)$ to
the log-joint. The posterior is derived [on its own page](gaussian.md).

### Class prior

In the original model the class probabilities $\pi_k$ are fixed. With `cell_type_prior`
set to `weighted`, $\pi$ is a random variable with a Dirichlet prior,

$$
p(\pi) = \mathrm{Dir}(\pi \mid \delta_1, \dots, \delta_K), \qquad
\log p(\pi) = \sum_k (\delta_k - 1) \log \pi_k + \text{const},
$$

over the real classes, with the Zero class held at a fixed weight. This adds
$\log p(\pi)$ to the log-joint. The posterior is derived [on its own page](class-prior.md).

### The full log-joint

The log-joint of the full model is the sum of equations (1) and (2) plus the priors on
the cells' geometry and on the class probabilities:

$$
\begin{aligned}
\log p(x, g, z, \zeta, \gamma, \eta, \theta, \rho, m, \Sigma, \pi) =
& - \sum_{g,c,k} \zeta_{c,k}\, \theta_c\, \mu_{g,k}\, A_c\, \gamma_{g,c}\, \eta_g \\
& + \sum_{s,c,k} z_{s,c}\, \zeta_{c,k}\,
   \log\!\big[\, \theta_c\, \mu_{g_s,k}\, e^{-D_c(x_s)}\, \gamma_{g_s,c}\, \eta_{g_s} \big] \\
& + \sum_{g,c} \log p(\gamma_{g,c}) + \sum_g \log p(\eta_g) + \sum_c \log p(\theta_c) \\
& + \sum_c \log p(m_c, \Sigma_c) \\
& + \sum_g \Big[ - \rho_g\, A_{\text{total}} + \sum_{s:\, g_s = g} z_{s,0}\, \log \rho_g
   + \log p(\rho_g) \Big] \\
& + \sum_{c,k} \zeta_{c,k}\, \log \pi_k + \sum_k (\delta_k - 1) \log \pi_k
   + \beta \sum_{c,k} \sum_{c' \in \mathcal{N}_c} \mathbf{1}(\zeta_{c,k} = \zeta_{c',k} = 1) .
\end{aligned}
\tag{3}
$$

## The variational approximation

The posterior is intractable and is approximated by **variational inference**: the member
of a tractable, factorised family closest to the true posterior in Kullback-Leibler
divergence, fitted by **coordinate ascent** (CAVI), one factor at a time. $\gamma$ depends
on $\theta$ and on the cell's class, so the three form one structured factor:

$$
p(z, \zeta, \gamma, \eta, \theta, \rho, m, \Sigma, \pi \mid x, g)
\approx q(\gamma \mid \zeta, \theta)\, q(\theta \mid \zeta)\, q(\zeta)\, q(z)\, q(\eta)\, q(\rho)\,
q(m, \Sigma)\, q(\pi) .
$$

Coordinate ascent updates each factor to its optimal form, which is always the **expected
log-joint over all the other factors**:

$$
\log q^*(x_j) = \mathbb{E}_{q(\text{rest})}\big[\log p(x, g, z, \zeta, \gamma, \eta, \theta, \rho, m, \Sigma, \pi)\big] + \text{const} .
$$

Each derivation that follows applies this equation: the terms of the
[log-joint](#the-generative-model) that involve the factor $x_j$ are kept, the expectation
over the other factors is taken, and the posterior is read off. $\theta_c$ is the
exception: as a point estimate it is **maximised** rather than integrated, against the same
expected log-joint.

The per-cell scale $\theta_c$ enters the intensity multiplicatively with $\gamma_{g,c}$.
With both as full random variables the marginalisation has no closed form and the Negative
Binomial likelihood of the cell typing is lost. $q(\theta_c)$ is therefore restricted to a
**point estimate** (a Dirac delta), which keeps $\theta_c$ constant during the update for
$\gamma_{g,c}$.

## Derivations

- **[Misread density $\rho_g$](misread-density.md)** - the per-gene background rate.
- **[The scale factors $\theta_c$, $\gamma_{g,c}$, $\eta_g$](scale-factors.md)** - the
  per-cell, per-gene-per-cell, and per-gene corrections to the reference.
- **[Cell-class assignment $q(\zeta)$](cell-class.md)** - the Negative Binomial likelihood,
  the class prior, and the MRF spatial term.
- **[Cell geometry $q(m_c, \Sigma_c)$](gaussian.md)** - the centre and
  covariance of each cell's Gaussian, and the expected log density the spot assignment
  uses.
- **[Spot-to-cell assignment $q(z)$](spot-assignment.md)** - which cell (or the background)
  each spot is assigned to.
- **[Class prior $q(\pi)$](class-prior.md)** - the Dirichlet posterior of the class
  probabilities under `weighted`.