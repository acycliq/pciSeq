---
id: cell-scale-theta
title: Derivation - the cell scale factor
sidebar_label: Cell scale factor
sidebar_position: 3
---

# Derivation: the cell scale factor $\theta_c$

The cell scale factor $\theta_c$ is an extension to the original model. It captures the
fact that some cells simply yield more transcripts than the reference predicts and others
fewer, applying a single whole-cell correction across all of a cell's genes. We give it a
conjugate prior $\theta_c \sim \mathrm{Gamma}(r_\theta, r_\theta)$, which centres the
correction at a baseline of $1.0$.

## Why a point estimate

$\theta_c$ enters the intensity multiplicatively with the per-gene-per-cell factor
$\gamma_{g,c}$. If both were treated as full random variables, marginalising the
fluctuations would require integrating over the product of two Gamma distributions, giving
a Poisson-Gamma-Gamma mixture with no closed form, and the Negative Binomial likelihood
that drives cell typing would break. We therefore restrict the variational distribution of
$\theta_c$ to a point estimate (a Dirac delta centred at $\theta_c^*$):

$$
q(\theta_c) = \delta(\theta_c - \theta_c^*) .
$$

The Dirac makes the expectation exact, $\mathbb{E}[\log\theta_c] = \log\mathbb{E}[\theta_c]$,
so $\theta_c$ is a constant during the update for $\gamma_{g,c}$ and the Poisson-Gamma
conjugacy is preserved. Optimising the point estimate against the log-joint *including* the
prior term makes this a **Maximum A Posteriori (MAP)** step: the prior regularises the
estimate, shrinking it toward the baseline $1.0$ when the data for a cell are sparse. The
status of mixing a Dirac factor with full variational factors is set out in the
[self-consistency appendix](appendix-self-consistency.md).

## The objective

This is a maximisation step. Conditioning on the class ($\zeta_{c,k} = 1$), the terms of
the log-joint that depend on $\theta_c$ are

$$
\mathcal{L}(\theta_c) =
- \theta_c \sum_g \mu_{g,k}\, A_c\, \gamma_{g,c}\, \eta_g
+ \sum_s z_{s,c}\, \log\theta_c
+ \underbrace{(r_\theta - 1)\log\theta_c - r_\theta\,\theta_c}_{\text{prior}}
+ \text{const}.
$$

Before differentiating we take the expectation over the other variational distributions
$q(z)$, $q(\gamma)$, and $q(\eta)$, replacing the latent variables by their expected
values:

$$
\mathbb{E}_{\gamma,\eta,z}[\mathcal{L}(\theta_c)]
= - \theta_c \Big( A_c \sum_g \mu_{g,k}\, \bar\gamma_{g,c}\, \bar\eta_g + r_\theta \Big)
  + \log\theta_c \Big( \sum_s \bar z_{s,c} + r_\theta - 1 \Big) + \text{const}.
$$

## Solving

Differentiating with respect to $\theta_c$ and setting the derivative to zero, with
$\bar N_c = \sum_s \bar z_{s,c}$ the expected total gene count of cell $c$:

$$
- \Big( A_c \sum_g \mu_{g,k}\, \bar\gamma_{g,c}\, \bar\eta_g + r_\theta \Big)
+ \frac{\bar N_c + r_\theta - 1}{\theta_c} = 0 .
$$

Solving gives the estimate

$$
\boxed{\;
\hat\theta_c = \frac{\bar N_c + r_\theta - 1}
                    {r_\theta + \sum_g A_c\, \mu_{g,k}\, \bar\gamma_{g,c}\, \bar\eta_g}
\;}
$$

## Reading the result

Up to the prior terms, the numerator is the observed number of spots in cell $c$, and the
denominator is the total count predicted for it under class $k$. So $\hat\theta_c \mid k$
is again an **observed-over-expected** ratio: greater than $1$ when the cell captures more
transcripts than the model predicts, less than $1$ when it captures fewer. A large
$r_\theta$ pulls the estimate toward the prior mean $1$, shrinking the correction when the
evidence is weak. Like $\gamma_{g,c}$, the estimate is computed conditional on the class
$k$, so it needs no weighting by $\bar\zeta_{c,k}$ in its own update.
