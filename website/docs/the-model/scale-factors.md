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
answers how efficiently gene $g$ is detected across the whole experiment. Its conjugate
Gamma posterior, with $N_g$ the total spots of gene $g$, is

$$
\boxed{\;
q(\eta_g)
= \mathrm{Gamma}\Big(
    r_\eta + N_g,\;\;
    \frac{r_\eta}{\eta_0}
    + \sum_{c,k} \bar\zeta_{c,k}\, \mu_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\theta_c
  \Big)
\;}
$$

Here the rate carries the baseline efficiency $\eta_0$ (typically $0.2$) through the prior
term $r_\eta/\eta_0$, and the sum over cells and classes is weighted by the soft class
assignments $\bar\zeta_{c,k}$, since $\eta_g$ is shared and must aggregate the evidence
from every cell. This is the **reparameterised** efficiency described in
[errata item 5](errata.md): the constant $\eta_0$ is pulled out of the prior and into the
intensity, leaving $\eta_g$ as a relative scaling factor with prior mean $1.0$.

## The three scale factors at a glance

| Factor | Indexed by | Corrects |
| --- | --- | --- |
| $\theta_c$ | cell (per candidate class) | the whole cell's total count |
| $\gamma_{g,c}$ | gene and cell (given the class) | one gene in one cell |
| $\eta_g$ | gene only (global) | one gene across the whole experiment |

The [warping the reference](../how-it-works/warping-the-reference.md) page gives the
intuition for how these stack from broad to fine.
