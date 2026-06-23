
# Derivation: the spot-to-cell assignment $q(z)$

The last block of the loop assigns each spot to the cell most likely to have produced it,
or to the background. The latent variable is the indicator $z_{s,c}$, which is $1$ when spot
$s$ belongs to cell $c$. We derive its variational posterior $q(z_{s,c})$ the same way as
the other factors: keep the terms of the log-joint that involve $z_{s,c}$, take expectations
over everything else, and read off the result.

## Which terms involve $z_{s,c}$

Spots are modelled as a spatial Poisson process with intensity
$\lambda_{g,c}(x) = \theta_c\, \mu_{g,k(c)}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g$. A Poisson
log-likelihood has two parts, $-\!\int\!\lambda(x)\,dx + \sum_s \log\lambda(x_s)$. The
indicator $z_{s,c}$ appears **only in the second part**: it picks out the log-intensity of
the cell each spot is assigned to. Keeping just those terms,

$$
\log p(\dots) \supset
\sum_{s,c,k} z_{s,c}\, \zeta_{c,k}\,
  \log\!\big[\theta_c\, \mu_{g_s,k}\, e^{-D_c(x_s)}\, \gamma_{g_s,c}\, \eta_{g_s}\big]
\;+\; \sum_s z_{s,0}\, \log\rho_{g_s} ,
$$

where the last sum is the background option ($c = 0$), whose intensity is the per-gene
[misread density](misread-density.md) $\rho_{g_s}$.

The first Poisson part, $-\!\int\!\lambda\,dx$ (the expected total count, and the equivalent
$\rho_g A_{\text{total}}$ for the background), contains **no** $z_{s,c}$. It is the same
whichever cell the spot is handed to, so it is constant across the assignment and drops out
under the normalisation below. This is why the area of the tissue never enters the
competition.

## The update for a cell ($c > 0$)

The coordinate-ascent update sets $\log q(z_{s,c}=1)$ to the expectation of those terms over
the other factors:

$$
\log q(z_{s,c}=1)
= \mathbb{E}_{\zeta,\theta,\gamma,\eta}\!\Big[
  \sum_k \zeta_{c,k}\,
  \big(\log\theta_c - D_c(x_s) + \log\mu_{g_s,k} + \log\gamma_{g_s,c} + \log\eta_{g_s}\big)
\Big] + \text{const}.
$$

Now carry the expectations inside, and mind the **class conditioning**. The cell scale
$\theta_c$ and the gene-cell factor $\gamma_{g_s,c}$ are both estimated *conditional on the
class* (see [scale-theta](scale-theta.md) and [scale-gamma](scale-gamma.md)), so inside the
$k$-th term they take their class-$k$ values:
$\mathbb{E}[\log\theta_c] = \log\bar\theta_{c\mid k}$ and
$\mathbb{E}[\log\gamma_{g_s,c}] = \overline{\log\gamma}_{g_s,c\mid k}$. The efficiency
$\mathbb{E}[\log\eta_{g_s}] = \overline{\log\eta}_{g_s}$ is gene-only, with no $k$; and
$\mathbb{E}[\zeta_{c,k}] = \bar\zeta_{c,k}$ is the cell-class posterior. The distance
$D_c(x_s)$ is constant and $\sum_k \bar\zeta_{c,k} = 1$, so it comes out of the sum:

$$
\log q(z_{s,c}=1)
= - D_c(x_s)
  + \sum_k \bar\zeta_{c,k}
    \big(\log\bar\theta_{c\mid k} + \overline{\log\gamma}_{g_s,c\mid k} + \overline{\log\eta}_{g_s} + \log\mu_{g_s,k}\big)
  + \text{const}.
$$

Exponentiating gives the unnormalised score $S_{s,c}$ for assigning the spot to cell $c$.

## The background option ($c = 0$)

The only term carrying $z_{s,0}$ is the background log-intensity $\log\rho_{g_s}$ (its area
term $\rho_g A_{\text{total}}$ is constant in $z$ and cancels with the rest). So

$$
\log q(z_{s,0}=1) = \overline{\log\rho_{g_s}} + \text{const},
$$

the **expected log** misread density of the spot's gene, $\overline{\log\rho_{g_s}} =
\psi(\hat r_{g_s}) - \log\hat\beta_{g_s}$ (the digamma form derived on the
[misread density](misread-density.md) page). Exponentiated, the background contributes
$\exp(\overline{\log\rho_{g_s}})$ to the competition, with no distance, scale, or efficiency
term attached. Note this is **not** the posterior mean rate $\bar\rho_{g_s} = \mathbb{E}[\rho_{g_s}]$:
because $\mathbb{E}[\log\rho] \neq \log\mathbb{E}[\rho]$, the term $\exp(\overline{\log\rho_{g_s}})$
sits strictly below $\bar\rho_{g_s}$. This is the exact quantity the code uses for the
background column (`genes.log_rho_bar`).

## Normalisation

The indicator $z_s$ picks exactly one option, so the scores are normalised across the nearby
cells and the background by a softmax:

$$
q\big(c(s)=c\big) = \frac{\exp(S_{s,c})}{\sum_{c'>0}\exp(S_{s,c'}) + \exp(\overline{\log\rho_{g_s}})},
\qquad
q\big(c(s)=0\big) = \frac{\exp(\overline{\log\rho_{g_s}})}{\sum_{c'>0}\exp(S_{s,c'}) + \exp(\overline{\log\rho_{g_s}})} .
$$

Dropping the shared normaliser, the posterior is

$$
\boxed{\;
q\big(c(s)=c\big)
\propto
\begin{cases}
\exp\!\Big[
  - D_c(x_s)
  + \displaystyle\sum_k
    \Big(
      \bar\zeta_{c,k} \log\bar\theta_{c\mid k}
      + \bar\zeta_{c,k} \overline{\log\gamma}_{g_s,c\mid k}
      + \bar\zeta_{c,k} \overline{\log\eta}_{g_s}
      + \bar\zeta_{c,k} \log\mu_{g_s,k}
    \Big)
\Big], & c > 0, \\[2pt]
\exp(\overline{\log\rho_{g_s}}), & c = 0 \ \text{(background)} .
\end{cases}
\;}
$$

## Reading the terms

The score for a cell $c > 0$ adds one **quantitative** term (geometry) to four
**qualitative** terms (expression); the background option $c = 0$ scores
$\exp(\overline{\log\rho_{g_s}})$. Each term is named and read below - just enough to see
what the symbol *means*. The conceptual narrative, the diagram, and the worked two-cell
examples live on the [how-it-works page](../how-it-works/spots-to-cells.md).

**Spatial fit** ($-D_c(x_s)$). The Gaussian log-likelihood of the spot's position under the
cell: $D_c(x)$ is the (Mahalanobis) distance to the cell's centre under its Gaussian shape,
so $e^{-D_c(x)}$ is the Gaussian weight and $-D_c(x_s)$ its logarithm. This is the only
quantitative term - a hard geometric measure of how well the spot sits inside the cell's
footprint, blind to which gene it carries. Nearer spots score higher.

**Alignment** ($\sum_k \bar\zeta_{c,k}\log\mu_{g_s,k}$). An inner product between the cell's
class posterior $\bar\zeta_c$ (how confident we are about its type) and the gene's expression
profile $\log\mu_{g_s,\cdot}$ (which types express the gene). Since $\bar\zeta_c$ sums to $1$
it equals $\mathbb{E}_{k\sim\bar\zeta_c}[\log\mu_{g_s,k}]$, the gene's expected log-expression
under the cell's own belief about its class. It is large only when both line up: a confident
type that also expresses the gene. *(class $\leftrightarrow$ gene)*

**Gravity** ($\sum_k \bar\zeta_{c,k}\log\bar\theta_{c\mid k}$). The same confidence-weighted
average, now of $\bar\theta_{c\mid k}$, the cell's total observed count over what class $k$
predicts. A cell capturing more transcripts than its type expects has $\bar\theta_{c\mid k} >
1$, a sparse one below $1$. Read it as the cell's mass: a heavier cell pulls harder, so all
else equal a spot drifts toward whichever cell is already capturing the most. *(cell size)*

**Enrichment** ($\sum_k \bar\zeta_{c,k}\log\bar\gamma_{g_s,c\mid k}$). The same form again,
now of $\bar\gamma_{g_s,c\mid k}$, *this cell's* observed-over-expected for the gene. Easily
confused with the alignment but distinct: alignment is class $\leftrightarrow$ gene (the
type's stereotype, shared by every cell of that type), enrichment is cell $\leftrightarrow$
gene (this individual cell's departure from its type). Because the rate factorises as
$\mu \times \gamma$, enrichment is exactly the residual the alignment leaves unexplained - it
is what tells two same-type cells apart. *(cell $\leftrightarrow$ gene)*

**Misread correction** ($\overline{\log\eta}_{g_s}$). The gene's detection efficiency.
Gene-only, so it is identical for every cell and cancels in any cell-versus-cell comparison;
it bites only against the **background**, attenuating a poorly detected gene's signal so its
spots are more readily called misreads (see the [note below](#the-efficiency-term-and-the-signal-to-noise-ratio)).

## The efficiency term and the signal-to-noise ratio

A subtle point, and the subject of [errata item 1](errata.md): although the efficiency
term $\overline{\log\eta}_{g_s}$ is the same for every cell $c > 0$, it does **not** cancel
during normalisation, because the assignment is also compared against the background
$\rho_{g_s}$, which carries no efficiency term. A low-efficiency gene therefore has its
signal attenuated relative to the background, making its spots more likely to be deemed
misreads. The original paper omitted this term, effectively treating every gene as perfectly
detected during assignment; including it lets $\eta$ act as a gene-specific scaling of the
signal-to-noise ratio.
