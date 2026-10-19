
# The scale factors

Three corrections sit between the cell type definitions and the observed counts, each
rescaling the expected expression at a different level of detail. The
[warping the cell type definitions](../how-it-works/warping-the-reference.md) page
describes them; this page derives them.

| Factor | Indexed by | Corrects |
| --- | --- | --- |
| **[$\theta_{c\mid k}$](#theta)** | cell, class | the whole cell's total count, given the class |
| **[$\gamma_{g,c\mid k}$](#gamma)** | gene, cell, class | a given gene in a given cell, under a given class |
| **[$\eta_g$](#eta)** | gene only (global) | one gene across the whole experiment |

$\theta_c$ is one number for a cell, $\gamma_{g,c}$ one per gene in that cell, and
$\eta_g$ one per gene shared by every cell. Each absorbs the mismatch at its own scale.

$\gamma_{g,c}$ is kept as a full Gamma random variable so it can be **integrated out**,
which collapses the Poisson count model into the Negative Binomial the
[cell-class assignment](cell-class.md) scores against. $\theta_c$ is held as a point
estimate for that reason: a second mixing distribution in the rate would not collapse.

## The derivations

### Derivation: the cell scale factor $\theta_c$ {#theta}

The cell scale factor $\theta_c$, also called the cell inefficiency, is an extension to
the original model. Some cells yield more transcripts than the reference predicts and
others fewer; $\theta_c$ is a single
whole-cell correction across all of a cell's genes. Its prior is the conjugate
$\theta_c \sim \mathrm{Gamma}(r_\theta, r_\theta)$, mean $1$.

**Conditional on the class.** Both $\theta_c$ and the cell-gene factor
[$\gamma_{g,c}$](#gamma) are computed **conditional on a candidate class $k$**. The
predicted count a cell is compared against depends on which class the cell is assumed to be,
so the correction is class-specific. The estimate is written $\hat\theta_{c\mid k}$, and the
[cell-class assignment](cell-class.md) recomputes it for every class $k$ it tests the cell
against.

#### The point estimate

$\theta_c$ enters the intensity multiplicatively with the per-gene-per-cell factor
$\gamma_{g,c}$. If both were treated as full random variables, marginalising the
fluctuations would require integrating over the product of two Gamma distributions, giving
a Poisson-Gamma-Gamma mixture with no closed form, and the Negative Binomial likelihood
that drives cell typing would be lost. The variational distribution of $\theta_c$ is
therefore restricted to a point estimate (a Dirac delta centred at $\theta_c^*$):

$$
q(\theta_c) = \delta(\theta_c - \theta_c^*) .
$$

The Dirac makes the expectation exact, $\mathbb{E}[\log\theta_c] = \log\mathbb{E}[\theta_c]$,
so $\theta_c$ is a constant during the update for $\gamma_{g,c}$ and the Poisson-Gamma
conjugacy is preserved. Optimising the point estimate against the log-joint *including* the
prior term makes this a **Maximum A Posteriori (MAP)** step: the prior regularises the
estimate, shrinking it toward the baseline $1.0$ when the data for a cell are sparse.

#### The objective

This is a maximisation step. Conditioning on the class ($\zeta_{c,k} = 1$), the terms of
the log-joint that depend on $\theta_c$ are

$$
\mathcal{L}(\theta_c) =
- \theta_c \sum_g \mu_{g,k}\, A_c\, \gamma_{g,c}\, \eta_g
+ \sum_s z_{s,c}\, \log\theta_c
+ \underbrace{(r_\theta - 1)\log\theta_c - r_\theta\,\theta_c}_{\text{prior}}
+ \text{const}.
$$

The expectation over the other variational distributions $q(z)$, $q(\gamma)$ and $q(\eta)$
replaces the latent variables by their expected values:

$$
\mathbb{E}_{\gamma,\eta,z}[\mathcal{L}(\theta_c)]
= - \theta_c \Big( A_c \sum_g \mu_{g,k}\, \bar\gamma_{g,c}\, \bar\eta_g + r_\theta \Big)
  + \log\theta_c \Big( \sum_s \bar z_{s,c} + r_\theta - 1 \Big) + \text{const}.
$$

#### Solving

Differentiating with respect to $\theta_c$ and setting the derivative to zero, with
$\bar N_c = \sum_s \bar z_{s,c}$ the expected total gene count of cell $c$:

$$
- \Big( A_c \sum_g \mu_{g,k}\, \bar\gamma_{g,c}\, \bar\eta_g + r_\theta \Big)
+ \frac{\bar N_c + r_\theta - 1}{\theta_c} = 0 .
$$

Solving gives the estimate

$$
\boxed{\;
\hat\theta_{c\mid k} = \frac{\bar N_c + r_\theta - 1}
                    {r_\theta + \sum_g A_c\, \mu_{g,k}\, \bar\gamma_{g,c}\, \bar\eta_g}
\;}
\tag{1}
$$

The class $k$ enters through the reference mean $\mu_{g,k}$ in the denominator: a different
candidate class gives a different predicted total, and hence a different $\hat\theta_{c\mid k}$.

#### Reading the result

Up to the prior terms, the numerator is the observed number of spots in cell $c$, and the
denominator is the total count predicted for it under class $k$. So $\hat\theta_{c\mid k}$
is again an **observed-over-expected** ratio: greater than $1$ when the cell captures more
transcripts than the model predicts, less than $1$ when it captures fewer. A large
$r_\theta$ moves the estimate toward the prior mean $1$, reducing the correction when the
evidence is weak. Like $\gamma_{g,c}$, the estimate is computed conditional on the class
$k$, so it needs no weighting by $\bar\zeta_{c,k}$ in its own update.

#### The prior strength $r_\theta$ {#r-theta}

$r_\theta$ sets how far the data can move $\theta_c$ from $1$. The two limits:

- **Weak prior ($r_\theta$ small).** $\hat\theta_{c\mid k} \approx \bar N_c / (\text{predicted
  total under } k)$, the cell's observed count over what the class predicts.
- **Strong prior ($r_\theta \to \infty$).** $\hat\theta_{c\mid k} \to 1$, and the
  whole-cell correction is switched off.

A weak prior has one consequence for **near-empty cells**. For a cell with few spots,
small $\bar N_c$, $\hat\theta_{c\mid k}$ can collapse to a small value, and a small
$\theta$ scales the predicted expression of any class $k$ down toward the cell's counts. A
real class shrunk far enough predicts almost nothing, the same as an empty cell.

The **Zero class**, which expects no expression, is the class for such cells (see
[assigning cells to classes](../how-it-works/cell-to-celltype.md)). When the prior is
weak and $\theta$ collapses, the cell can instead be assigned to a real class scaled down
to predict almost nothing, and the cell gets a spurious class.

Why the shrunk class wins comes out of the estimate. The total count the model predicts for
cell $c$ under class $k$ is $\theta$ times the sum over genes of the expected counts,

$$
\theta_c \sum_g A_c\, \mu_{g,k}\, \bar\gamma_{g,c}\, \bar\eta_g \;=\; \theta_c\, E_{c\mid k},
$$

writing $E_{c\mid k}$ for the sum, the same sum as in the denominator of equation (1).
Replace $\theta_c$ by the estimate of equation (1):

$$
\hat\theta_{c\mid k}\, E_{c\mid k}
\;=\; (\bar N_c + r_\theta - 1)\,\frac{E_{c\mid k}}{E_{c\mid k} + r_\theta}.
$$

When $r_\theta$ is small against $E_{c\mid k}$ the last fraction is close to $1$, and the
predicted total is about $\bar N_c + r_\theta - 1$ whatever the class: every class predicts
the cell's own count plus a little. The total no longer separates the classes; what remains
is how each class splits it over the genes, so the few spots decide the class by which genes
they are, not by how many there are.

An example at $r_\theta = 2$. A cell holds $0.48$ spots, one Plp1 spot at probability
$0.44$. Before scaling, Oligo expects $356$ spots in the cell and Zero expects $55$.

| | $\hat\theta$ at $r_\theta = 2$ | prediction | $\hat\theta$ at $r_\theta = 25$ | prediction |
|---|---|---|---|---|
| Oligo | $0.0041$ | $1.47$ | $0.064$ | $22.9$ |
| Zero | $0.026$ | $1.43$ | $0.306$ | $16.8$ |

At $r_\theta = 2$ both classes predict about $1.5$ spots, Oligo expects that spot to be Plp1,
and Oligo wins at $0.94$. At $r_\theta = 25$ Oligo predicts $23$ spots for half a spot, loses
on the amount, and the cell goes to Zero.

The same $\theta$ enters the [spot assignment](spot-assignment.md#reading-the-terms):
$\sum_k \bar\zeta_{c,k}\log\hat\theta_{c\mid k}$ is added to the score of every spot the cell
competes for. The background has no such term. For the cell above, as Oligo at
$r_\theta = 2$, the term is $\log 0.0041 \approx -5.5$, and the one spot it holds sits at
probability $0.44$, close to level with the background. At $r_\theta = 25$ the cell is Zero
and the term is $\log 0.306 \approx -1.2$.

The weak prior has a use too. A real cell with poor capture, a thin slice or a cell at the
edge of the section, holds far fewer spots than its class usually gives, and a free $\theta$
still finds its class from the composition of those spots.

A weaker prior lets the data set $\theta$ at the risk of near-empty cells being explained
by a collapsed class; a stronger prior holds $\theta$ near $1$ so those cells go to Zero.
$r_\theta$ acts as a pseudo-count against the cell's own spots: at $2$ a cell with ten spots
sets its own scale, at $25$ it needs of the order of $25$. The default `rTheta` is $25$.

### Derivation: the cell-gene scale factor $\gamma_{g,c}$ {#gamma}

The cell-gene scale factor $\gamma_{g,c}$ adjusts a single gene in a single cell. It is the
most fine-grained of the three [scale factors](#the-scale-factors), but it also plays a
special structural role: it is the factor that turns the Poisson count model into a
**Negative Binomial**, and that Negative Binomial is what the
[cell-class assignment](cell-class.md) scores each cell against. Unlike the cell scale
$\theta_c$, which is a [point estimate](#theta), $\gamma_{g,c}$ is kept as a full
Gamma random variable and **integrated out**.

**Conditional on the class.** Like [$\theta_c$](#theta), $\gamma_{g,c}$ is computed
**conditional on a candidate class $k$**: it measures how far gene $g$ in cell $c$ deviates
from the expression that class $k$ predicts, so the answer depends on which class the cell is
assumed to be. Every formula on this page carries that conditioning, written
$\gamma_{g,c}\mid k$, and the [cell-class assignment](cell-class.md) recomputes it for each
class it tests.

#### The count is a Poisson-Gamma mixture

Fix a gene $g$ and a cell $c$ of class $k$. The expected number of spots of that gene in
that cell, before the per-gene-cell fluctuation, is the deterministic predicted count

$$
\lambda_{g,c} = \mu_{g,k}\, A_c\, \bar\eta_g\, \bar\theta_c ,
$$

the reference expression $\mu_{g,k}$ scaled by the inside-cell-bonus factor $A_c$ (see the
[notation table](overview.md#notation); $A_c = 1$ when the bonus is off), the gene efficiency
$\bar\eta_g$ and the cell scale $\bar\theta_c$. The scale factor $\gamma_{g,c}$ multiplies
this rate, and the observed count is Poisson around it:

$$
N_{c,g} \mid \gamma_{g,c} \;\sim\; \mathrm{Poisson}\big(\lambda_{g,c}\, \gamma_{g,c}\big),
\qquad
\gamma_{g,c} \sim \mathrm{Gamma}(r_\gamma, r_\gamma) .
$$

The prior $\mathrm{Gamma}(r_\gamma, r_\gamma)$ has mean $1$, so a priori the count sits at
$\lambda_{g,c}$ and deviates cell by cell through $\gamma_{g,c}$. A Poisson whose rate is a
Gamma random variable is a **Poisson-Gamma mixture**, which accommodates the overdispersion
of transcript counts.

#### Deriving the posterior $q(\gamma_{g,c})$

By the mean-field (CAVI) update, the optimal factor for $\gamma_{g,c}$ is the expected
log-joint over every other latent, with the $\gamma_{g,c}$-free terms folded into the
constant:

$$
\log q^*(\gamma_{g,c}) = \mathbb{E}_{-\gamma_{g,c}}\big[\log p(x,g,z,\zeta,\gamma,\eta,\theta)\big] + \text{const}.
$$

Conditioning on the class assignments $\zeta$, the class prior, the spatial term, and
$\bar N_{c,g}\,\overline{\log\theta_c}$ are all constant in $\gamma_{g,c}$. Three pieces
survive: the spatial integral contributes the linear term $-\lambda_{g,c}\,\gamma_{g,c}$, the
$\bar N_{c,g}$ expected spots each contribute a $\log\gamma_{g,c}$, and the prior
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

#### Integrating $\gamma$ out: the Negative Binomial

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
variance exceeds the mean, which a Poisson cannot produce; the Gamma mixing is what
accommodates the **overdispersion** of real counts. The shape $r_\gamma$ sets how much:
small $r_\gamma$ allows large deviations, large $r_\gamma$ approaches a Poisson at
$\lambda_{g,c}$. How large is relative to the expected count: the excess variance
$\lambda_{g,c}^2/r_\gamma$ is small next to $\lambda_{g,c}$ only when
$r_\gamma \gg \lambda_{g,c}$, so a gene the class barely expects is close to a Poisson
already at $r_\gamma = 2$, while a strongly expressed gene remains overdispersed at
values far above it.

This Negative Binomial is the per-gene likelihood that the
[cell-class assignment](cell-class.md) multiplies across genes to score a cell against each
candidate class. Keeping $\gamma_{g,c}$ conjugate, so that it can be integrated out in closed
form, is therefore what keeps cell typing tractable. It is also why $\theta_c$ is held as a
[point estimate](#theta): if it were a full random variable too, the rate would be a product
of two mixing distributions and this clean Poisson-Gamma collapse would be lost.

### Derivation: the in situ efficiency $\eta_g$ {#eta}

The efficiency $\eta_g$ is the most global of the three [scale factors](#the-scale-factors): a
per-gene parameter shared across all cells and classes, answering how efficiently gene $g$
is detected across the whole experiment. It is estimated as a **relative** factor $\eta_g'$
with prior mean $1.0$, while the baseline detection rate $\eta_0$ (the `Inefficiency`
setting, default $0.2$) stays as an explicit constant next to the reference mean - the
intensity carries the product $\eta_0\,\mu_{g,k}$, not a renamed symbol:

$$
\lambda_{g,c}(x) = \theta_c\,\eta_0\,\mu_{g,k(c)}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g' .
$$

#### Per gene, shared across cells

Detection efficiency differs between genes: probe chemistry, sequence and length make some
transcripts easier to read out than others, so $\eta_g$ is estimated **per gene**. It is
shared across **all** cells because the detection rate of a gene is a property of the
assay, not of a cell. $\gamma_{g,c}$ varies cell by cell: $\eta_g$ is how well the gene is
detected anywhere, $\gamma_{g,c}$ how far this gene in this cell deviates from what its class and
$\theta_c$ predict.

#### Deriving the posterior $q(\eta_g')$

The prior is $\eta_g' \sim \mathrm{Gamma}(r_\eta, r_\eta)$, whose log-density is

$$
\log p(\eta_g') = (r_\eta - 1)\log\eta_g' - r_\eta\,\eta_g' .
$$

By the mean-field (CAVI) update, the optimal factor for $\eta_g'$ is the expected log-joint
over every other latent, with all the $\eta_g'$-free terms folded into the constant:

$$
\log q^*(\eta_g') = \mathbb{E}_{-\eta_g'}\big[\log p(x,g,z,\zeta,\gamma,\eta,\theta)\big] + \text{const}.
$$

Writing $\eta_g = \eta_0\,\eta_g'$, only three pieces of the joint carry $\eta_g'$. The
predicted-count term is linear in $\eta_g'$; each observed spot of gene $g$ contributes a
$\log\eta_g'$ through the per-spot log-intensity (since $\log\eta_g = \log\eta_0 + \log\eta_g'$
and $\log\eta_0$ is constant); and the prior adds its two terms:

$$
\log q(\eta_g')
= -\Big(\sum_{c,k}\zeta_{c,k}\,\theta_c\,\eta_0\mu_{g,k}\,A_c\,\gamma_{g,c}\Big)\eta_g'
  + \Big(\sum_{s:\,g_s=g}\sum_{c,k} z_{s,c}\,\zeta_{c,k}\Big)\log\eta_g'
  + (r_\eta-1)\log\eta_g' - r_\eta\,\eta_g'
  + \text{const}.
$$

Now simplify the $\log\eta_g'$ coefficient, the double sum over the gene-$g$ spots and over
the cell-class pairs. Both $z_{s,c}$ and $\zeta_{c,k}$ are latent indicators: a cell has one
true class and a spot one true parent cell, but neither is known, so inference carries the
posterior probabilities over them, the soft assignments $\bar z_{s,c}$ and the class
probabilities $\bar\zeta_{c,k}$. Taking the expectation splits the product under the mean-field
factorisation, $\mathbb{E}[z_{s,c}\,\zeta_{c,k}] = \bar z_{s,c}\,\bar\zeta_{c,k}$, and the
sum over $k$ uses that a cell's class probabilities form a distribution,
$\sum_k \bar\zeta_{c,k} = 1$. The coefficient becomes

$$
\bar N_g = \sum_{s:\,g_s=g}\sum_c \bar z_{s,c},
$$

the **expected** number of reads of gene $g$ assigned to cells. Replacing the remaining latents
by their variational means ($\gamma\!\to\!\bar\gamma$, $\theta\!\to\!\bar\theta$) and gathering
the $\log\eta_g'$ and the linear $\eta_g'$ terms:

$$
\log q(\eta_g')
= (\bar N_g + r_\eta - 1)\log\eta_g'
  - \Big(r_\eta + \sum_{c,\,k\neq\text{zero}} \bar\zeta_{c,k}\, \eta_0\mu_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\theta_c\Big)\eta_g'
  + \text{const}.
$$

This is the log of a Gamma density; reading off its shape and rate:

$$
\boxed{\;
q(\eta_g')
= \mathrm{Gamma}\Big(
    \bar N_g + r_\eta,\;\;
    r_\eta + \sum_{c,\,k\neq\text{zero}} \bar\zeta_{c,k}\, \eta_0\mu_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\theta_c
  \Big)
\;}
$$

*Note.* Background and Zero-cell reads are not produced by a class profile, so in practice $\bar N_g$
excludes them. $\eta_g'$ is then reads observed in real-class cells (ie any class except Zero) over reads predicted for them.

The sum runs over all cells and candidate classes, weighted by the soft class assignments
$\bar\zeta_{c,k}$, since $\eta_g'$ is shared and aggregates evidence from every cell. A value
$\eta_g' > 1$ means gene $g$ is detected better than the baseline, $\eta_g' < 1$ worse.

#### The efficiency term against the background

$\eta_g$ depends only on the gene, so it takes the same value for every candidate cell in
a [spot-to-cell assignment](spot-assignment.md) and cancels between cells. It acts in one
place: the comparison between a cell and the **background**. The background scores a spot
with the [misread density](misread-density.md) $\rho_g$ alone, which carries no efficiency
term, so $\eta_g$ is the one factor on the cell's side without a counterpart. A gene with a
small $\eta_g$ has its signal attenuated against the background and its spots are more
readily called misreads.

The original paper dropped this term, which inflated the signal-to-noise ratio for poorly
detected genes.
