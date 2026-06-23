
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

### The typo

As printed, the paper's efficiency prior is $\eta_g \sim \mathrm{Gamma}(r_\eta, \eta_0)$ with
$\eta_0 = 0.2$ ($r_\eta$ is the shape, the `rGene` setting, default $20$). In the shape-rate
parameterisation this has mean $r_\eta/\eta_0 = 100$, not $0.2$. The intended prior, with
mean $\eta_0$, is $\eta_g \sim \mathrm{Gamma}(r_\eta, r_\eta/\eta_0)$.

### The reparameterisation

The whole construction rests on the **scale property of the Gamma distribution**:

$$
X \sim \mathrm{Gamma}(a, \beta)
\quad\Longrightarrow\quad
cX \sim \mathrm{Gamma}(a, \beta/c).
$$

Rather than estimate the absolute efficiency, the implementation keeps the baseline $\eta_0$
(the `Inefficiency` setting passed to `pciSeq.fit()`) as an explicit constant next to the
reference mean, $\eta_0\,\mu_{g,k}$, and estimates a **relative** factor
$\eta_g' = \eta_g/\eta_0$ with prior mean $1.0$. Collecting the $\eta_g'$ terms of the
expected log-joint (full derivation on the [efficiency](scale-factors.md#eta) page) gives

$$
\eta_g' \sim \mathrm{Gamma}\big(N_g + r_\eta,\; r_\eta + S_g\big),
\qquad
S_g = \sum_{c,k} \bar\zeta_{c,k}\, \eta_0\mu_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\theta_c .
$$

The relative factor reads off cleanly against the baseline: $\eta_g' = 1$ is a gene detected
at exactly the assumed rate, $\eta_g' > 1$ better than the baseline, $\eta_g' < 1$ worse.
Pulling $\eta_0$ out of the prior also avoids the typo above: the prior on $\eta_g'$ is
simply $\mathrm{Gamma}(r_\eta, r_\eta)$, whose mean is $1$ by construction, so there is no
baseline constant left in the prior to get wrong.

The scale property links the two parameterisations: with $c = \eta_0$, the absolute variable
$\eta_g = \eta_0\,\eta_g'$ recovers the prior $\mathrm{Gamma}(r_\eta, r_\eta/\eta_0)$ and
divides the posterior rate by $\eta_0$. The two are the same model written two ways; the
relative one is just better conditioned for the optimiser.

### Summary

| | Prior | Posterior |
| --- | --- | --- |
| Relative, $\eta_g'$ | $\mathrm{Gamma}(r_\eta, r_\eta)$ | $\mathrm{Gamma}\big(N_g + r_\eta,\ r_\eta + S_g\big)$ |
| Absolute, $\eta_g = \eta_0\,\eta_g'$ | $\mathrm{Gamma}(r_\eta, r_\eta/\eta_0)$ | $\mathrm{Gamma}\big(N_g + r_\eta,\ (r_\eta + S_g)/\eta_0\big)$ |

with $N_g$ the observed spots of gene $g$ and $S_g$ as above (the per-cell factor
$\bar\theta_c$ is this model's extension; the original paper omits it). The code uses the
relative form: it centres the working variable at $1.0$, preconditions the optimisation near
the right order of magnitude, and vectorises cleanly across genes.
