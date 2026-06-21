---
id: scale-factors
title: Derivation - the scale factors gamma and eta
sidebar_label: Scale factors
sidebar_position: 4
---

# Derivation: the scale factors $\gamma_{g,c}$ and $\eta_g$

Two further corrections sit alongside the [cell scale factor](cell-scale-theta.md)
$\theta_c$. The per-gene-per-cell factor $\gamma_{g,c}$ adjusts a single gene in a single
cell; the per-gene efficiency $\eta_g$ is a global correction shared across all cells. Both
keep their conjugate Gamma form.

## The cell-gene scale factor $\gamma_{g,c}$

Conditioning on the class assignments $\zeta$, the class prior and the spatial term are
constant in $\gamma$, and so is $\bar N_{c,g}\,\overline{\log\theta_c}$. Collecting the
remaining terms of the joint and adding the Gamma prior
$\log p(\gamma_{g,c}) = (r_\gamma - 1)\log\gamma_{g,c} - r_\gamma\,\gamma_{g,c}$:

$$
\log q(\gamma_{g,c} \mid k(c))
= (\bar N_{c,g} + r_\gamma - 1)\log\gamma_{g,c}
  - \big(r_\gamma + \mu_{g,k(c)}\, A_c\, \bar\eta_g\, \bar\theta_c\big)\gamma_{g,c}
  + \text{const}.
$$

This is the log of a Gamma density with shape $\bar N_{c,g} + r_\gamma$ and rate
$r_\gamma + \mu_{g,k(c)}\, A_c\, \bar\eta_g\, \bar\theta_c$. Therefore

$$
\boxed{\;
q(\gamma_{g,c} \mid \theta, k(c))
= \mathrm{Gamma}\big(
    \bar N_{c,g} + r_\gamma,\;\;
    r_\gamma + \mu_{g,k(c)}\, A_c\, \bar\eta_g\, \bar\theta_c
  \big)
\;}
$$

The posterior mean

$$
\bar\gamma_{g,c}
= \frac{\bar N_{c,g} + r_\gamma}{r_\gamma + \mu_{g,k}\, A_c\, \bar\eta_g\, \bar\theta_c}
$$

has the same observed-over-expected structure as $\hat\theta_c$, but for a *single* gene:
the numerator is the expected count of gene $g$ in cell $c$, the denominator the prediction
for that same gene. Where $\theta_c$ corrects the whole cell at once, $\gamma_{g,c}$
corrects each gene individually.

## The in situ efficiency $\eta_g$

The efficiency $\eta_g$ is a per-gene parameter shared across all cells and classes: it
answers how efficiently gene $g$ is detected across the whole experiment. It is estimated as
a **relative** factor $\eta_g'$ with prior mean $1.0$, while the baseline detection rate
$\eta_0$ (the `Inefficiency` setting, default $0.2$) stays as an explicit constant next to
the reference mean - the intensity carries the product $\eta_0\,\mu_{g,k}$, not a renamed
symbol:

$$
\lambda_{g,c}(x) = \eta_0\,\mu_{g,k(c)}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g' .
$$

The prior is $\eta_g' \sim \mathrm{Gamma}(r_\eta, r_\eta)$, whose log-density is

$$
\log p(\eta_g') = (r_\eta - 1)\log\eta_g' - r_\eta\,\eta_g' .
$$

Keep the terms of the expected log-joint that involve $\eta_g'$. The spatial integral
contributes the predicted-count sum (a linear term in $\eta_g'$), the $N_g$ observed spots of
gene $g$ each contribute a $\log\eta_g'$, and the prior adds the two terms above:

$$
\log q(\eta_g')
= N_g \log\eta_g'
  - \Big(\sum_{c,k} \bar\zeta_{c,k}\, \eta_0\mu_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\theta_c\Big)\eta_g'
  + \underbrace{(r_\eta - 1)\log\eta_g' - r_\eta\,\eta_g'}_{\text{prior } \log p(\eta_g')}
  + \text{const}.
$$

Now gather the $\log\eta_g'$ terms and the linear $\eta_g'$ terms separately:

$$
\log q(\eta_g')
= (N_g + r_\eta - 1)\log\eta_g'
  - \Big(r_\eta + \sum_{c,k} \bar\zeta_{c,k}\, \eta_0\mu_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\theta_c\Big)\eta_g'
  + \text{const}.
$$

This is the log of a Gamma density; reading off its shape and rate:

$$
\boxed{\;
q(\eta_g')
= \mathrm{Gamma}\Big(
    N_g + r_\eta,\;\;
    r_\eta + \sum_{c,k} \bar\zeta_{c,k}\, \eta_0\mu_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\theta_c
  \Big)
\;}
$$

The sum runs over all cells and candidate classes, weighted by the soft class assignments
$\bar\zeta_{c,k}$, since $\eta_g'$ is shared and aggregates evidence from every cell. A value
$\eta_g' > 1$ means gene $g$ is detected better than the baseline, $\eta_g' < 1$ worse. The
equivalent absolute form $\eta_g = \eta_0\,\eta_g'$, and the prior/posterior summary, are in
[errata item 5](errata.md).

## The three scale factors at a glance

| Factor | Indexed by | Corrects |
| --- | --- | --- |
| $\theta_c$ | cell (per candidate class) | the whole cell's total count |
| $\gamma_{g,c}$ | gene and cell (given the class) | one gene in one cell |
| $\eta_g$ | gene only (global) | one gene across the whole experiment |

The [warping the reference](../how-it-works/warping-the-reference.md) page gives the
intuition for how these stack from broad to fine.
