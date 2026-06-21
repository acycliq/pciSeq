---
id: scale-gamma
title: Derivation - the cell-gene scale factor gamma
sidebar_label: gamma - cell-gene scale
sidebar_position: 3
---

# Derivation: the cell-gene scale factor $\gamma_{g,c}$

The cell-gene scale factor $\gamma_{g,c}$ adjusts a single gene in a single cell. It is the
most fine-grained of the three [scale factors](scale-factors.md), but it also plays a
special structural role: it is the factor that turns the Poisson count model into a
**Negative Binomial**, and that Negative Binomial is what the
[cell-class assignment](cell-class.md) scores each cell against. Unlike the cell scale
$\theta_c$, which is a [point estimate](scale-theta.md), $\gamma_{g,c}$ is kept as a full
Gamma random variable and **integrated out**.

**Conditional on the class.** Like [$\theta_c$](scale-theta.md), $\gamma_{g,c}$ is computed
**conditional on a candidate class $k$**: it measures how far gene $g$ in cell $c$ deviates
from the expression that class $k$ predicts, so the answer depends on which type the cell is
assumed to be. Every formula on this page carries that conditioning, written
$\gamma_{g,c}\mid k$, and the [cell-class assignment](cell-class.md) recomputes it for each
class it tests.

## The count is a Poisson-Gamma mixture

Fix a gene $g$ and a cell $c$ of class $k$. The expected number of spots of that gene in
that cell, before the per-gene-cell fluctuation, is the deterministic predicted count

$$
\lambda_{g,c} = \mu_{g,k}\, A_c\, \bar\eta_g\, \bar\theta_c ,
$$

the reference expression $\mu_{g,k}$ scaled by the cell area $A_c$, the gene efficiency
$\bar\eta_g$ and the cell scale $\bar\theta_c$. The scale factor $\gamma_{g,c}$ multiplies
this rate, and the observed count is Poisson around it:

$$
N_{c,g} \mid \gamma_{g,c} \;\sim\; \mathrm{Poisson}\big(\lambda_{g,c}\, \gamma_{g,c}\big),
\qquad
\gamma_{g,c} \sim \mathrm{Gamma}(r_\gamma, r_\gamma) .
$$

The prior $\mathrm{Gamma}(r_\gamma, r_\gamma)$ has mean $1$, so a priori the count sits at
$\lambda_{g,c}$; the gene is free to deviate cell by cell through $\gamma_{g,c}$. A Poisson
whose own rate is a Gamma random variable is a **Poisson-Gamma mixture**, and this is the
structure that makes the model robust to the overdispersion real transcript counts show.

## Deriving the posterior $q(\gamma_{g,c})$

Conditioning on the class assignments $\zeta$, the class prior, the spatial term, and
$\bar N_{c,g}\,\overline{\log\theta_c}$ are all constant in $\gamma_{g,c}$. Keep the terms of
the expected log-joint that involve it: the spatial integral contributes the linear term
$-\lambda_{g,c}\,\gamma_{g,c}$, the $\bar N_{c,g}$ spots each contribute a
$\log\gamma_{g,c}$, and the prior
$\log p(\gamma_{g,c}) = (r_\gamma - 1)\log\gamma_{g,c} - r_\gamma\,\gamma_{g,c}$ adds two more:

$$
\log q(\gamma_{g,c} \mid k)
= \bar N_{c,g}\log\gamma_{g,c}
  - \lambda_{g,c}\,\gamma_{g,c}
  + \underbrace{(r_\gamma - 1)\log\gamma_{g,c} - r_\gamma\,\gamma_{g,c}}_{\text{prior } \log p(\gamma_{g,c})}
  + \text{const}.
$$

Gather the $\log\gamma_{g,c}$ terms and the linear $\gamma_{g,c}$ terms separately:

$$
\log q(\gamma_{g,c} \mid k)
= (\bar N_{c,g} + r_\gamma - 1)\log\gamma_{g,c}
  - (r_\gamma + \lambda_{g,c})\,\gamma_{g,c}
  + \text{const}.
$$

This is the log of a Gamma density; reading off its shape and rate:

$$
\boxed{\;
q(\gamma_{g,c} \mid k(c))
= \mathrm{Gamma}\big(
    \bar N_{c,g} + r_\gamma,\;\;
    r_\gamma + \mu_{g,k}\, A_c\, \bar\eta_g\, \bar\theta_c
  \big)
\;}
$$

The posterior is conjugate, as the Poisson-Gamma structure guarantees. Its mean

$$
\bar\gamma_{g,c}
= \frac{\bar N_{c,g} + r_\gamma}{r_\gamma + \mu_{g,k}\, A_c\, \bar\eta_g\, \bar\theta_c}
$$

is the **observed-over-expected** ratio for a single gene: the expected count
$\bar N_{c,g}$ over the prediction $\lambda_{g,c}$, regularised by $r_\gamma$. Where
$\theta_c$ rescales the whole cell at once, $\gamma_{g,c}$ rescales each gene individually.
Both are computed conditional on the class $k$.

## Integrating $\gamma$ out: the Negative Binomial

The reason $\gamma_{g,c}$ is kept as a full Gamma rather than a point estimate is what
happens when it is **marginalised**. For the Poisson-Gamma pair above, integrating the rate
$\gamma_{g,c}$ against its $\mathrm{Gamma}(r_\gamma, r_\gamma)$ prior gives a closed-form
Negative Binomial for the count:

$$
N_{c,g} \mid \gamma_{g,c} \sim \mathrm{Poisson}(\lambda_{g,c}\,\gamma_{g,c}),
\quad
\gamma_{g,c} \sim \mathrm{Gamma}(r_\gamma, r_\gamma)
\;\;\Longrightarrow\;\;
N_{c,g} \sim \mathrm{NB}(r_\gamma, \lambda_{g,c}),
$$

with mean $\lambda_{g,c}$ and variance $\lambda_{g,c} + \lambda_{g,c}^2/r_\gamma$. The
variance exceeds the mean, which a plain Poisson could never produce: the Gamma mixing is
exactly what lets the model accommodate the **overdispersion** of real counts. The shape
$r_\gamma$ controls how much: small $r_\gamma$ allows large deviations, large $r_\gamma$
pulls the count back toward a pure Poisson at $\lambda_{g,c}$.

This Negative Binomial is the per-gene likelihood that the
[cell-class assignment](cell-class.md) multiplies across genes to score a cell against each
candidate type. Keeping $\gamma_{g,c}$ conjugate, so that it can be integrated out in closed
form, is therefore what keeps cell typing tractable. It is also why $\theta_c$ and the
latent shift $\mathbf{b}_c$ are held as [point estimates](scale-theta.md): if they were full
random variables too, the rate would be a product of several mixing distributions and this
clean Poisson-Gamma collapse would be lost.
