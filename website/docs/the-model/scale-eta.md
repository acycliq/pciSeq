---
id: scale-eta
title: Derivation - the in situ efficiency eta
sidebar_label: eta - gene efficiency
sidebar_position: 4
---

# Derivation: the in situ efficiency $\eta_g$

The efficiency $\eta_g$ is the most global of the three [scale factors](scale-factors.md): a
per-gene parameter shared across all cells and classes, answering how efficiently gene $g$
is detected across the whole experiment. It is estimated as a **relative** factor $\eta_g'$
with prior mean $1.0$, while the baseline detection rate $\eta_0$ (the `Inefficiency`
setting, default $0.2$) stays as an explicit constant next to the reference mean - the
intensity carries the product $\eta_0\,\mu_{g,k}$, not a renamed symbol:

$$
\lambda_{g,c}(x) = \eta_0\,\mu_{g,k(c)}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g' .
$$

## Why per-gene, and why global

Detection efficiency is not uniform across genes: probe chemistry, sequence, and length
make some transcripts far easier to read out than others. A single global efficiency would
force the noisy and the clean genes to share one number, so $\eta_g$ is estimated **per
gene**. But it is shared across **all** cells, because the detection rate of a gene is a
property of the assay, not of any one cell. That is what separates it from $\gamma_{g,c}$,
which varies cell by cell: $\eta_g$ asks "how well is this gene read out anywhere?", while
$\gamma_{g,c}$ asks "how far does this gene in this cell deviate from its class?".

## Deriving the posterior $q(\eta_g')$

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

## What $\eta_g$ decides: cell versus background

Because $\eta_g$ depends only on the gene, it takes the same value for every candidate cell
in a [spot-to-cell assignment](spot-assignment.md). When two genuine cells compete for a
spot, $\eta_g$ is a common offset on both sides and cancels: it has no say in which cell
wins. Its influence shows up in exactly one place - the contest between a cell and the
**background**.

The background claims a spot through a single quantity: the
[misread density](misread-density.md) $\rho_g$, a score that carries no efficiency term at
all. So $\eta_g$ is the one factor on the cell's side with no counterpart on the
background's side, and that is precisely why it survives the comparison. It sets, gene by
gene, how strong a cell's signal must be to win a spot rather than have it written off as
noise. A low-efficiency gene has a small $\eta_g$, so its signal is attenuated against the
background and its spots are more readily called misreads - as they should be, since that
gene really is detected poorly.

The original paper dropped this term, which inflated the signal-to-noise ratio for poorly
detected genes; the correction is [errata item 1](errata.md).
