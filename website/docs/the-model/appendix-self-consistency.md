
# Appendix: self-consistency of the variational construction

This appendix records why the construction is sound: how the updates behave under a Dirac
variational family, the standard identities they rely on, and what objective is actually
being optimised. Nothing here changes the model or the algorithm.

## The optimum under a Dirac variational family

Mean-field variational inference approximates the intractable posterior by the member of a
tractable family closest in Kullback-Leibler divergence. Our family is **structured**: the
scale factors $\gamma_{g,c}$ and $\eta_g$ and the class assignments $\zeta_{c,k}$ keep full
conjugate forms (Gamma, Gamma, and Categorical), while the per-cell scale $\theta_c$ is
restricted to a Dirac point estimate:

$$
q(\theta_c \mid k) = \delta(\theta_c - \hat\theta_{c\mid k}) .
$$

Under a Dirac, every expectation is exact by the sifting property:
$\mathbb{E}_q[f(\theta_c)\mid k] = f(\hat\theta_{c\mid k})$ for any measurable $f$. In
particular the expectation of the log equals the log of the expectation,

$$
\mathbb{E}_q[\theta_c \mid k] = \hat\theta_{c\mid k},
\qquad
\mathbb{E}_q[\log\theta_c \mid k] = \log\hat\theta_{c\mid k} = \log\mathbb{E}_q[\theta_c\mid k] ,
$$

which is exactly what lets $\theta_c$ pass through the $\gamma$ update as a constant and keep
the Poisson-Gamma conjugacy intact. These are exact properties of the Dirac measure, not
approximations. No Jensen-type inequality applies, because Jensen requires a non-degenerate
measure and the Dirac is degenerate.

## Standard identities used

- **Poisson-Gamma as Negative Binomial.** If $N \mid \gamma \sim \mathrm{Poisson}(K\gamma)$
  and $\gamma \sim \mathrm{Gamma}(r, r)$, then marginally $N \sim \mathrm{NB}(r, K)$ with
  mean $K$ and variance $K + K^2/r$ - the collapse that drives the cell-class update.
- **Gamma log-moment.** If $X \sim \mathrm{Gamma}(a, b)$ then
  $\mathbb{E}[\log X] = \psi(a) - \log b$, with $\psi$ the digamma function - the source of
  every $\overline{\log\cdot}$ term in the updates ($\overline{\log\gamma}$,
  $\overline{\log\eta}$, $\overline{\log\rho}$).

## The objective: Variational EM, not pure VB

The Evidence Lower Bound contains an entropy term $-\mathbb{E}_q[\log q]$. For a Dirac
factor this differential entropy is not finite, so the ELBO is not literally well-defined
when the family contains a point mass.

The resolution is standard. With a Dirac component the scheme is **Variational EM**: the
Dirac variable is treated as a parameter and point-estimated by maximising the expected
log-joint, while the non-Dirac variables keep their full variational posteriors. The
objective is the free energy

$$
\mathcal{F}\big(\{\hat{\boldsymbol\beta}_j\}, \{q_i\}\big)
= \mathbb{E}_{\prod_i q_i}\big[\log p(\text{data, latents}, \{\hat{\boldsymbol\beta}_j\})\big]
  + \sum_i H[q_i],
$$

where the point estimate is $\hat\theta_{c\mid k}$ and the $\{q_i\}$ are the non-Dirac
factors. The point estimate contributes no entropy term because it is a parameter, not a
distribution. Each conjugate update and each MAP step performs coordinate ascent on
$\mathcal{F}$, and the model treats $\theta_c$ in exactly this way.
