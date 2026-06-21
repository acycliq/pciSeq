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

The key point is **where** this term actually changes anything. The efficiency
$\overline{\log\eta}_{g_s}$ depends only on the spot's gene, not on the cell, so it takes
the **same value for every candidate cell** $c > 0$. In a comparison between two genuine
neighbouring cells it is a common offset that cancels: it never affects which cell wins.
It matters only in the **spot-versus-background** decision. The background option ($c = 0$)
carries no efficiency term at all, so $\overline{\log\eta}_{g_s}$ is exactly the asymmetry
between "assign to some cell" and "assign to the background". Dropping it, as the original
expression does, removes that asymmetry and lets a low-efficiency gene compete against the
background as if it were perfectly detected, which is precisely the signal-to-noise
distortion above.

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

### The original form

As written in the paper, the gene efficiency has the prior

$$
\eta_g \sim \mathrm{Gamma}(r, \eta_0),
\qquad \eta_0 = 0.2 .
$$

This is a typo. In the rate parameterisation it has mean $r/\eta_0 = 100$ (for $r = 20$),
not $0.2$. To give the intended mean $\eta_0 = 0.2$ the prior should be

$$
\eta_g \sim \mathrm{Gamma}(r, r/\eta_0),
\qquad \mathbb{E}[\eta_g] = \frac{r}{r/\eta_0} = \eta_0 = 0.2 .
$$

So as printed the prior does not have the intended mean.

### The reparameterisation

In the paper's (absolute) parameterisation the variational posterior for the efficiency is
$q(\eta_g) = \mathrm{Gamma}(r_\eta + N_g,\ r_\eta/\eta_0 + S_g)$, so its posterior mean is

$$
\mathbb{E}[\eta_g]
= \frac{r_\eta + N_g}{\dfrac{r_\eta}{\eta_0} + S_g},
$$

with $N_g$ the observed spots of gene $g$ and $S_g$ the rate sum defined in the table below.

The reparameterisation rests on the **scale property of the Gamma distribution**: for
$c > 0$,

$$
X \sim \mathrm{Gamma}(a, \beta)
\quad\Longrightarrow\quad
cX \sim \mathrm{Gamma}(a, \beta/c),
$$

that is, scaling the variable leaves the shape untouched and divides the rate by $c$. With
$c = \eta_0$ it lets us move between the absolute efficiency $\eta_g$ and the rescaled
working variable $\eta_g'$ (defined below):

$$
\eta_g' \sim \mathrm{Gamma}(r_\eta, r_\eta)
\quad\Longrightarrow\quad
\eta_g = \eta_0\,\eta_g' \sim \mathrm{Gamma}\!\big(r_\eta,\ r_\eta/\eta_0\big),
$$

which is exactly the absolute prior from the original form above.

The implementation does not patch the prior in place. Instead it **pulls the baseline
constant $\eta_0$ out of the prior**, floats it as a free factor in the intensity function,
and **fuses it into the reference mean** $\mu_{g,k}$. Starting from the intensity with the
baseline written out explicitly,

$$
\lambda_{g,c}(x) =
\underbrace{\eta_0 \cdot \mu_{g,k(c)}}_{\text{adjusted expression } \mu'_{g,k}}
\cdot\; e^{-D_c(x)} \cdot \gamma_{g,c} \cdot \eta_g' .
$$

The product $\eta_0\,\mu_{g,k}$ is the **adjusted (pre-scaled) expression**
$\mu'_{g,k} = \eta_0\,\mu_{g,k}$: the scRNA-seq reference mean already discounted by the
baseline detection rate. Here $\eta_0$ is the **`Inefficiency`** setting the user passes in
the options dictionary to `pciSeq.fit()` (default $0.2$, i.e. a 20% baseline detection
rate). With the baseline absorbed into $\mu'$, the estimated variable $\eta_g'$ is no longer
an absolute efficiency but a **relative** scaling factor, centred at a prior mean of $1.0$:

$$
\eta_g' \sim \mathrm{Gamma}(r_\eta, r_\eta), \qquad \mathbb{E}[\eta_g'] = 1 .
$$

A value $\eta_g' > 1$ means gene $g$ is detected better than the 20% baseline, $\eta_g' < 1$
worse. With $\eta_0$ baked into $\mu'$ this way, the paper's posterior equations hold as
written (subject to item 1), reading $\mu_{g,k}$ as $\mu'_{g,k}$ and $\eta_g$ as the
relative factor $\eta_g'$.

### Side by side

| | Without reparameterisation | With reparameterisation |
| --- | --- | --- |
| Estimated variable | $\eta_g$ - absolute efficiency, prior mean $\eta_0 \approx 0.2$ | $\eta_g'$ - relative factor, prior mean $1.0$ |
| Reference mean | $\mu_{g,k}$ (raw scRNA-seq) | $\mu'_{g,k} = \eta_0\,\mu_{g,k}$ (pre-scaled) |
| Prior on the variable | $\mathrm{Gamma}(r_\eta,\ r_\eta/\eta_0)$ | $\mathrm{Gamma}(r_\eta,\ r_\eta)$ |
| Intensity $\lambda_{g,c}(x)$ | $\mu_{g,k}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g$ | $\mu'_{g,k}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g'$ |
| Posterior $q(\eta)$ | $\mathrm{Gamma}\!\big(r_\eta + N_g,\ \tfrac{r_\eta}{\eta_0} + S_g\big)$ | $\mathrm{Gamma}\!\big(r_\eta + N_g,\ r_\eta + S_g'\big)$ |

where $N_g$ is the total observed spots of gene $g$ and the rate sum runs over all cells and
candidate classes,

$$
S_g = \sum_{c,k} \bar\zeta_{c,k}\, \mu_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\theta_c,
\qquad
S_g' = \sum_{c,k} \bar\zeta_{c,k}\, \mu'_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\theta_c
     = \eta_0\, S_g .
$$

The two columns describe the **same model**: substituting $\eta_g = \eta_0\,\eta_g'$ and
$\mu'_{g,k} = \eta_0\,\mu_{g,k}$ turns one into the other. Only the bookkeeping differs -
the right-hand column keeps the working variable centred at $1.0$, which is what makes it
better behaved numerically. The same scaling carries the posterior across too, since
$S_g' = \eta_0\,S_g$ makes the right-hand rate $r_\eta + S_g'$ equal to $\eta_0$ times the
left-hand rate $r_\eta/\eta_0 + S_g$.

This is numerically superior: it preconditions the optimisation near the right order of
magnitude (the 20% baseline), gives the shape parameter $r_\eta$ an intuitive reading as
pseudo-observations of that baseline, and vectorises cleanly across genes.
