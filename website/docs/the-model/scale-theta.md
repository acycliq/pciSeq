
# Derivation: the cell scale factor $\theta_c$

The cell scale factor $\theta_c$ is an extension to the original model. It captures the
fact that some cells simply yield more transcripts than the reference predicts and others
fewer, applying a single whole-cell correction across all of a cell's genes. We give it a
conjugate prior $\theta_c \sim \mathrm{Gamma}(r_\theta, r_\theta)$, which centres the
correction at a baseline of $1.0$.

**Conditional on the class.** Both $\theta_c$ and the cell-gene factor
[$\gamma_{g,c}$](scale-gamma.md) are computed **conditional on a candidate class $k$**. The
predicted count a cell is compared against depends on which type the cell is assumed to be,
so the correction is class-specific. We write the estimate $\hat\theta_{c\mid k}$, and the
[cell-class assignment](cell-class.md) recomputes it for every class $k$ it tests the cell
against.

## Why a point estimate

$\theta_c$ enters the intensity multiplicatively with the per-gene-per-cell factor
$\gamma_{g,c}$. If both were treated as full random variables, marginalising the
fluctuations would require integrating over the product of two Gamma distributions, giving
a Poisson-Gamma-Gamma mixture with no closed form, and the Negative Binomial likelihood
that drives cell typing would break. We therefore restrict the variational distribution of
$\theta_c$ to a point estimate (a Dirac delta centred at $\theta_c^*$):

$$
q(\theta_c) = \delta(\theta_c - \theta_c^*) .
$$

The Dirac makes the expectation exact, $\mathbb{E}[\log\theta_c] = \log\mathbb{E}[\theta_c]$,
so $\theta_c$ is a constant during the update for $\gamma_{g,c}$ and the Poisson-Gamma
conjugacy is preserved. Optimising the point estimate against the log-joint *including* the
prior term makes this a **Maximum A Posteriori (MAP)** step: the prior regularises the
estimate, shrinking it toward the baseline $1.0$ when the data for a cell are sparse. The
status of mixing a Dirac factor with full variational factors is set out in the
[self-consistency appendix](appendix-self-consistency.md).

## The objective

This is a maximisation step. Conditioning on the class ($\zeta_{c,k} = 1$), the terms of
the log-joint that depend on $\theta_c$ are

$$
\mathcal{L}(\theta_c) =
- \theta_c \sum_g \mu_{g,k}\, A_c\, \gamma_{g,c}\, \eta_g
+ \sum_s z_{s,c}\, \log\theta_c
+ \underbrace{(r_\theta - 1)\log\theta_c - r_\theta\,\theta_c}_{\text{prior}}
+ \text{const}.
$$

Before differentiating we take the expectation over the other variational distributions
$q(z)$, $q(\gamma)$, and $q(\eta)$, replacing the latent variables by their expected
values:

$$
\mathbb{E}_{\gamma,\eta,z}[\mathcal{L}(\theta_c)]
= - \theta_c \Big( A_c \sum_g \mu_{g,k}\, \bar\gamma_{g,c}\, \bar\eta_g + r_\theta \Big)
  + \log\theta_c \Big( \sum_s \bar z_{s,c} + r_\theta - 1 \Big) + \text{const}.
$$

## Solving

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
$$

The class $k$ enters through the reference mean $\mu_{g,k}$ in the denominator: a different
candidate type gives a different predicted total, and hence a different $\hat\theta_{c\mid k}$.

## Reading the result

Up to the prior terms, the numerator is the observed number of spots in cell $c$, and the
denominator is the total count predicted for it under class $k$. So $\hat\theta_{c\mid k}$
is again an **observed-over-expected** ratio: greater than $1$ when the cell captures more
transcripts than the model predicts, less than $1$ when it captures fewer. A large
$r_\theta$ pulls the estimate toward the prior mean $1$, shrinking the correction when the
evidence is weak. Like $\gamma_{g,c}$, the estimate is computed conditional on the class
$k$, so it needs no weighting by $\bar\zeta_{c,k}$ in its own update.

## The prior strength $r_\theta$

$r_\theta$ sets how far the data are allowed to move $\theta_c$ from its baseline of $1$.
Both extremes are legitimate; the right choice depends on how much you trust the counts:

- **Weak prior ($r_\theta$ small).** The posterior is essentially data-driven:
  $\hat\theta_{c\mid k} \approx \bar N_c / (\text{predicted total under } k)$, the cell's
  observed count over what the class predicts. The correction follows the data freely.
- **Strong prior ($r_\theta \to \infty$).** The posterior collapses onto the prior:
  $\hat\theta_{c\mid k} \to 1$, and the whole-cell correction is effectively switched off.

So a weak prior is a perfectly reasonable choice when you want the data to drive inference.
It just carries one pitfall worth keeping in mind, for **near-empty cells**.

Take a cell with very few spots (a small $\bar N_c$) under a weak prior. $\hat\theta_{c\mid k}$
is then free to collapse to a very small value, and a small $\theta$ scales the predicted
expression of **any** class $k$ down toward the cell's handful of counts. Shrink a real cell
type far enough and it "predicts" almost nothing, which is exactly what an empty cell looks
like.

The model already has a class for empty cells: the **Zero class**, which expects no
expression and absorbs debris and poorly segmented fragments (see
[assigning cells to cell types](../how-it-works/cell-to-celltype.md)). A near-empty cell
ought to land there. But when the prior is weak and $\theta$ collapses, a genuine type can be
shrunk down to imitate the Zero class and win the cell instead: the Zero class is skipped and
the cell gets a spurious type.

This is the tension to weigh when setting $r_\theta$. A weaker prior lets the data speak but
risks near-empty cells being explained away by a collapsed type; a stronger prior holds
$\theta$ near $1$ so those cells fall to the Zero class. The default is $25$ (the `rTheta`
setting), which leans toward the safe side - but it is a choice, not a rule.
