
# Derivation: the spot-to-cell assignment $q(z)$

The last step of the sweep assigns each spot to a cell or to the background. The latent
variable is the indicator $z_{s,c}$, which is $1$ when spot $s$ belongs to cell $c$. Its
variational posterior $q(z_{s,c})$ is derived as for the other factors: the terms of the
log-joint that involve $z_{s,c}$ are kept, the expectation over the rest is taken, and the
result is read off.

## The terms in $z_{s,c}$

Spots are modelled as a spatial Poisson process with intensity

$$
\lambda_{g,c}(x) = \theta_c\, \mu_{g,k(c)}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g .
$$

A Poisson log-likelihood has two parts:

$$
-\!\int\!\lambda(x)\,dx \;+\; \sum_s \log\lambda(x_s) .
$$

The indicator $z_{s,c}$ appears **only in the second part**: it picks out the log-intensity
of the cell each spot is assigned to. Keeping just those terms of the
[full log-joint](overview.md#the-generative-model),

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
whichever cell the spot is assigned to, so it is constant across the assignment and drops
out under the normalisation below. The extent of the tissue therefore does not enter the
assignment.

## The update for a cell ($c > 0$)

The [coordinate-ascent update](overview.md#the-variational-approximation) sets
$\log q(z_{s,c}=1)$ to the expectation of those terms over the other factors:

$$
\begin{aligned}
\log q(z_{s,c}=1)
&= \mathbb{E}_{\zeta,\theta,\gamma,\eta,m,\Sigma}\Big[ \sum_k \zeta_{c,k}\big(
   \log\theta_c - D_c(x_s) + \log\mu_{g_s,k} \\
&\qquad\qquad\quad {}+ \log\gamma_{g_s,c} + \log\eta_{g_s} \big) \Big] + \text{const}.
\end{aligned}
$$

Carrying the expectations inside, with the **class conditioning**: the cell scale
$\theta_c$ and the gene-cell factor $\gamma_{g_s,c}$ are both estimated *conditional on the
class* (see [scale-theta](scale-factors.md#theta) and [scale-gamma](scale-factors.md#gamma)), so inside the
$k$-th term they take their class-$k$ values:
$\mathbb{E}[\log\theta_c] = \log\bar\theta_{c\mid k}$ and
$\mathbb{E}[\log\gamma_{g_s,c}] = \overline{\log\gamma}_{g_s,c\mid k}$. The efficiency
$\mathbb{E}[\log\eta_{g_s}] = \overline{\log\eta}_{g_s}$ is gene-only, with no $k$; and
$\mathbb{E}[\zeta_{c,k}] = \bar\zeta_{c,k}$ is the cell-class posterior. The position term
is averaged over the cell's geometry, $\mathbb{E}_{m,\Sigma}[-D_c(x_s)] =
-\overline{D}_c(x_s)$, derived on the [cell geometry](gaussian.md) page. It
carries no $k$, and $\sum_k \bar\zeta_{c,k} = 1$, so the position and efficiency terms sit
outside the sum. Writing the assignment-dependent part as the score $S_{s,c}$, one named
term per line:

$$
\begin{aligned}
S_{s,c} = {}
& -\overline{D}_c(x_s) && \text{(position term)} \\
& + \sum_k \bar\zeta_{c,k}\,\log\mu_{g_s,k} && \text{(class expression)} \\
& + \sum_k \bar\zeta_{c,k}\,\log\bar\theta_{c\mid k} && \text{(cell scale)} \\
& + \sum_k \bar\zeta_{c,k}\,\overline{\log\gamma}_{g_s,c\mid k} && \text{(cell-gene scale)} \\
& + \overline{\log\eta}_{g_s} && \text{(gene efficiency)}
\end{aligned}
$$

so $\log q(z_{s,c}=1) = S_{s,c} + \text{const}$, and exponentiating, cell $c$ contributes
$\exp(S_{s,c})$ to the competition.

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
q\big(c(s)=c\big) = \frac{\exp(S_{s,c})}{Z},
\qquad
q\big(c(s)=0\big) = \frac{\exp(\overline{\log\rho_{g_s}})}{Z},
$$

with the shared normaliser $Z = \sum_{c'>0}\exp(S_{s,c'}) + \exp(\overline{\log\rho_{g_s}})$.
Dropping it, the posterior is

$$
\boxed{\;
q\big(c(s)=c\big)
\propto
\begin{cases}
\exp(S_{s,c}), & c > 0, \\[2pt]
\exp(\overline{\log\rho_{g_s}}), & c = 0 \ \text{(background)} .
\end{cases}
\;}
$$

## Reading the terms

The score for a cell $c > 0$ is one position term plus four expression terms. The
background option $c = 0$ scores $\overline{\log\rho_{g_s}}$. The
[how it works page](../how-it-works/spots-to-cells.md) gives the same terms without the
symbols.

**Position term** ($-\overline{D}_c(x_s)$). The expected log density of the spot's position
under the cell's Gaussian, normalisation included: minus half the squared Mahalanobis
distance to the posterior centre $\hat m_c$ with precision $\mathbb{E}[\Lambda_c]$, minus
the expected log normaliser, minus a flat cost $d / (2\hat\kappa_c)$ for the uncertainty
in the centre. The normaliser is what makes this term commensurate with the background
density $\rho_g$.
It does not depend on the gene. The closed form is on the
[cell geometry](gaussian.md#the-expected-log-density) page.

**Class expression** ($\sum_k \bar\zeta_{c,k}\log\mu_{g_s,k}$). The gene's log
expected expression averaged over the cell's class posterior $\bar\zeta_c$. It is large
when the cell is confidently of a type that expresses the gene. It is the same for every
cell of a given type.

**Cell scale** ($\sum_k \bar\zeta_{c,k}\log\bar\theta_{c\mid k}$). The same average of
$\bar\theta_{c\mid k}$, the cell's total observed count over what class $k$ predicts.
Above 1 for a cell holding more transcripts than its type expects, below 1 for a sparse
one. It does not depend on the gene.

**Cell-gene scale** ($\sum_k \bar\zeta_{c,k}\overline{\log\gamma}_{g_s,c\mid k}$). The
same average of the cell's observed-over-expected for this gene. Since the rate
factorises as $\mu \times \gamma$, this is the residual the class expression leaves, and
it is what separates two cells of the same type.

**Gene efficiency** ($\overline{\log\eta}_{g_s}$). The gene's detection efficiency.
It depends on the gene only, so it cancels between cells and acts only against the
background, which has no efficiency term. A poorly detected gene scores lower against the
background and its spots are more readily called misreads, see
[below](#the-efficiency-term-and-the-signal-to-noise-ratio).

## The efficiency term and the signal-to-noise ratio

Although the efficiency term $\overline{\log\eta}_{g_s}$ is the same for every cell
$c > 0$, it does **not** cancel
during normalisation, because the assignment is also compared against the background
$\rho_{g_s}$, which carries no efficiency term. A low-efficiency gene therefore has its
signal attenuated relative to the background, making its spots more likely to be deemed
misreads. The original paper omitted this term, effectively treating every gene as perfectly
detected during assignment; including it lets $\eta$ act as a gene-specific scaling of the
signal-to-noise ratio.
