---
id: spot-assignment
title: Derivation - the spot-to-cell assignment
sidebar_label: Spot-to-cell assignment
sidebar_position: 6
---

# Derivation: the spot-to-cell assignment $q(z)$

The last block of the loop assigns each spot to the cell most likely to have produced it,
or to the background. With the [cell types](cell-class.md) and the scale factors in hand,
the posterior over the assignment of spot $s$ is

$$
\boxed{\;
q\big(c(s)=c\big)
\propto
\begin{cases}
\exp\!\Big[
  - D_c(x_s)
  + \displaystyle\sum_k \bar\zeta_{c,k}
    \Big(
      \log\bar\theta_c
      + \overline{\log\gamma}_{g_s,c}
      + \overline{\log\eta}_{g_s}
      + \log\mu_{g_s,k}
    \Big)
\Big], & c > 0, \\[2pt]
\rho_g, & c = 0 \ \text{(background)} .
\end{cases}
\;}
$$

## Reading the terms

For a candidate cell $c > 0$ the score combines, in log space:

- $-D_c(x_s)$ - the **spatial term**: spots near the cell score higher;
- $\log\bar\theta_c$ - the cell's [overall scale](cell-scale-theta.md) (using
  $\overline{\log\theta_c} = \log\bar\theta_c$, exact under the point estimate);
- $\overline{\log\gamma}_{g_s,c}$ and $\overline{\log\eta}_{g_s}$ - the
  [gene-cell and gene efficiency](scale-factors.md) corrections;
- $\log\mu_{g_s,k}$ - the **compatibility** of the spot's gene with the cell's likely
  class, weighted by the cell-class posterior $\bar\zeta_{c,k}$.

For the background option $c = 0$ the score is the per-gene
[misread density](misread-density.md) $\rho_g$. A spot is assigned to the background unless
some nearby cell explains it better, which is how genuine misreads are filtered out.

## The efficiency term and the signal-to-noise ratio

A subtle point, and the subject of [errata item 1](errata.md): although the efficiency
term $\overline{\log\eta}_{g_s}$ is the same for every cell $c > 0$, it does **not** cancel
during normalisation, because the assignment is also compared against the background
$\rho_g$, which carries no efficiency term. A low-efficiency gene therefore has its signal
attenuated relative to the background, making its spots more likely to be deemed misreads.
The original paper omitted this term, effectively treating every gene as perfectly detected
during assignment; including it lets $\eta$ act as a gene-specific scaling of the
signal-to-noise ratio.
