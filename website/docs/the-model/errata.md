---
id: errata
title: Errata - corrections to Qian et al. (2020)
sidebar_label: Errata
sidebar_position: 8
---

# Errata: corrections to Qian et al. (2020)

This page records mathematical inconsistencies in the original pciSeq model
(Qian et al., 2020) and how the current model resolves them. They divide into functional
errors, typos in the text, and deliberate implementation strategies.

## 1. The missing efficiency term in spot assignment

The original variational update for the spot assignment $q(c(s)=c)$ reads

$$
q(c(s)=c) \propto \exp\Big[
  - D_c(x_s) + \overline{\log\gamma}_{g_s,c} + \sum_k \bar\zeta_{c,k}\log\mu_{g_s,k}
\Big].
$$

The term $\overline{\log\eta}_{g_s}$, the expected log detection efficiency, is **missing**.
It appears in the model's own intensity function, so taking the expectation of the
log-joint leaves it in the score for any cell $c > 0$.

Because the assignment is normalised against the background density $\rho$ (the $c=0$
class), which has no efficiency term, omitting $\eta$ does not cancel - it distorts the
signal-to-noise ratio. A low-efficiency gene ought to have its signal attenuated and be
more readily assigned to the background; the original expression instead treats every gene
as perfectly detected during assignment. The corrected expression restores the term:

$$
q(c(s)=c) \propto \exp\Big[
  - D_c(x_s) + \overline{\log\gamma}_{g_s,c} + \overline{\log\eta}_{g_s} + \sum_k \bar\zeta_{c,k}\log\mu_{g_s,k}
\Big].
$$

This is the form used on the [spot-to-cell assignment](spot-assignment.md) page.

## 2. Typo in the gamma prior

The text states $\gamma_{g,c} \sim \mathrm{Gamma}(r, 1)$ with shape $r = 2$. That gives
prior mean $\mathbb{E}[\gamma] = 2$, implying the model expects twice the scRNA-seq
expression by default. But the posterior rate update is $r + \mu_{g,k} A_c \bar\eta_g$,
which implies a prior rate of $r$, not $1$. The prior was intended to be
$\mathrm{Gamma}(r, r)$ (mean $1$), so the reference means $\mu$ are not incorrectly scaled.

## 3. Structured dependency in the efficiency update

The update for the gene efficiency $\eta$ uses a single expectation $\bar\gamma_{g,c}$. But
in the structured approximation $q(\zeta, \gamma) = q(\zeta)\,q(\gamma\mid\zeta)$ the
posterior for $\gamma$ is **conditional on the class** $k$, since $\gamma$ is a ratio of
observed to expected reads and the expectation depends on which class mean $\mu_{g,k}$ is
used. Using an unconditioned average biases the efficiency estimate.

## 4. Notation collision on $r$

The symbol $r$ denotes two independent quantities: the mean radius of the DAPI region
($r_{\text{DAPI}}$, page 8) and the dispersion parameter of the Negative Binomial
($r_{\text{NB}}$, pages 8-9).

## 5. The efficiency reparameterisation

The text states $\eta_g \sim \mathrm{Gamma}(r, \eta_0)$ with expected efficiency
$\eta_0 = 0.2$. In rate parameterisation that gives mean $r/\eta_0 = 100$, not $0.2$; for a
mean of $0.2$ the parameters should be $\mathrm{Gamma}(r, r/\eta_0)$.

The implementation resolves this by **reparameterisation**: it pulls the baseline constant
$\eta_0$ out of the prior and into the intensity function,

$$
\lambda_{g,c}(x) = \eta_0 \cdot \mu_{g,k(c)} \cdot e^{-D_c(x)} \cdot \gamma_{g,c} \cdot \eta_g',
$$

so the estimated variable $\eta_g'$ becomes a **relative** scaling factor with prior mean
$1.0$. This is numerically superior: it preconditions the optimisation near the right order
of magnitude (20%), gives the shape parameter an intuitive reading as pseudo-observations
of the baseline, and vectorises cleanly across genes. Equivalently, one can absorb $\eta_0$
into a pre-scaled reference mean $\mu'_{g,k} = \eta_0\,\mu_{g,k}$, after which the paper's
posterior equations hold as written (subject to item 1), with $\eta_g$ read as the relative
factor $\eta_g'$ and its prior set to $\mathrm{Gamma}(r_\eta, r_\eta)$.
