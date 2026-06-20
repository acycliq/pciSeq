---
id: appendix-self-consistency
title: Appendix - self-consistency of the construction
sidebar_label: Appendix - self-consistency
sidebar_position: 9
---

# Appendix: self-consistency of the variational construction

This appendix records why the construction is sound: how the updates behave under a Dirac
variational family, the standard identities they rely on, and what objective is actually
being optimised. Nothing here changes the model or the algorithm.

## The optimum under a Dirac variational family

Mean-field variational inference approximates the intractable posterior by the member of a
tractable family closest in Kullback-Leibler divergence. Our family is **structured**: the
scale factors $\gamma_{g,c}$ and $\eta_g$, the class covariance $\Sigma_k$, and the class
assignments $\zeta_{c,k}$ keep full conjugate forms (Gamma, Gamma, Inverse-Wishart, and
Categorical), while the per-cell scale $\theta_c$ and the latent shift $\mathbf{b}_c$ are
restricted to Dirac point estimates:

$$
q(\theta_c \mid k) = \delta(\theta_c - \hat\theta_{c\mid k}),
\qquad
q(\mathbf{b}_c \mid k) = \delta(\mathbf{b}_c - \hat{\mathbf{b}}_{c\mid k}) .
$$

Under a Dirac, every expectation is exact by the sifting property:
$\mathbb{E}_q[f(\mathbf{b}_c)\mid k] = f(\hat{\mathbf{b}}_{c\mid k})$ for any measurable
$f$. In particular

$$
\mathbb{E}_q[\mathbf{b}_c \mid k] = \hat{\mathbf{b}}_{c\mid k},
\qquad
\mathbb{E}_q[\mathbf{b}_c\mathbf{b}_c^{\mathrm{T}} \mid k] = \hat{\mathbf{b}}_{c\mid k}\hat{\mathbf{b}}_{c\mid k}^{\mathrm{T}},
\qquad
\mathbb{E}_q[e^{b_{g,c}} \mid k] = e^{\hat b_{g,c\mid k}} .
$$

These are exact properties of the Dirac measure, not approximations. No Jensen-type
inequality applies, because Jensen requires a non-degenerate measure and the Dirac is
degenerate.

## Standard identities used

- **Wishart / Inverse-Wishart duality.** If $\Sigma \sim \mathcal{IW}(\nu, \Psi)$ then
  $\Sigma^{-1} \sim \mathcal{W}(\nu, \Psi^{-1})$, so
  $\mathbb{E}[\Sigma^{-1}] = \nu\,\Psi^{-1}$ - the precision used in the Newton step for
  $\mathbf{b}_c$.
- **Inverse-Wishart log-determinant moment.** For dimension $p$,
  $\mathbb{E}[\log|\Sigma|] = \log|\Psi| - p\log 2 - \sum_{i=1}^p \psi\!\big(\tfrac{\nu - i + 1}{2}\big)$,
  with $\psi$ the digamma function - used inside the cell-class density term.
- **Inverse-Wishart conjugacy.** An $\mathcal{IW}(\nu_0, \Psi_0)$ prior with
  $\prod_i \mathcal{N}(\mathbf{b}_i \mid \mathbf{0}, \Sigma)$ gives posterior
  $\mathcal{IW}(\nu_0 + N, \Psi_0 + S)$ with scatter
  $S = \sum_i \mathbf{b}_i\mathbf{b}_i^{\mathrm{T}}$.
- **Poisson-Gamma as Negative Binomial.** If $N \mid \gamma \sim \mathrm{Poisson}(K\gamma)$
  and $\gamma \sim \mathrm{Gamma}(r, r)$, then marginally $N \sim \mathrm{NB}(r, K)$ with
  mean $K$ and variance $K + K^2/r$ - the collapse that drives the cell-class update.

## The objective: Variational EM, not pure VB

The Evidence Lower Bound contains an entropy term $-\mathbb{E}_q[\log q]$. For a Dirac
factor this differential entropy is not finite, so the ELBO is not literally well-defined
when the family contains a point mass.

The resolution is standard. With a Dirac component the scheme is **Variational EM**: the
Dirac variables are treated as parameters and point-estimated by maximising the expected
log-joint, while the non-Dirac variables keep their full variational posteriors. The
objective is the free energy

$$
\mathcal{F}\big(\{\hat{\boldsymbol\beta}_j\}, \{q_i\}\big)
= \mathbb{E}_{\prod_i q_i}\big[\log p(\text{data, latents}, \{\hat{\boldsymbol\beta}_j\})\big]
  + \sum_i H[q_i],
$$

where the $\{\hat{\boldsymbol\beta}_j\}$ are the point estimates ($\hat\theta_{c\mid k}$ and
$\hat{\mathbf{b}}_{c\mid k}$) and the $\{q_i\}$ are the non-Dirac factors. The point
estimates contribute no entropy term because they are parameters, not distributions. Each
conjugate update and each MAP step performs coordinate ascent on $\mathcal{F}$. The
original model treats $\theta_c$ in exactly this way; the gene-gene extension applies the
same construction unchanged to $\mathbf{b}_c$.

## Two expectations not to be conflated

The **prior-predictive mean** of a count, integrating $\gamma_{g,c}$ and $b_{g,c}$ against
their priors, carries a factor from the generative process:

$$
\mathbb{E}_{\text{prior-predictive}}[N_{c,g}] = K_{c,g}\, e^{(\Sigma_k)_{gg}/2},
\qquad K_{c,g} = \mu_{g,k}\, A_c\, \eta_g\, \theta_c .
$$

The **variational posterior expectation of the rate**, used inside every update, is a
different quantity under a different law:

$$
\mathbb{E}_q\big[K_{c,g}\, \gamma_{g,c}\, e^{b_{g,c}} \mid k\big]
= K_{c,g}\, \bar\gamma_{g,c\mid k}\, e^{\hat b_{g,c\mid k}} .
$$

There is no reason for the two to coincide, and the inference does not claim to compute the
first. The factor $e^{(\Sigma_k)_{gg}/2}$ is a moment of the generative model that any
inference scheme would face; it is not a bias of this one.
