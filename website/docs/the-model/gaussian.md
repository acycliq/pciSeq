
# Derivation: cell geometry $q(m_c, \Sigma_c)$

Each cell claims spots through a position term, the density of the spot's position under a
Gaussian centred on the cell. The geometry of a cell, on this page, is that Gaussian: its
centre $m_c$ and its covariance $\Sigma_c$, not the segmentation mask. In the original model
the Gaussian is fixed: its centre is the segmentation centroid and its covariance is
$r^2 I$, with $r$ the mean cell radius. Here $m_c$ and $\Sigma_c$ are random variables, one
pair per cell, with a prior and a variational posterior like every other factor. Spots that a cell is confident of
pull its centre towards them and stretch its covariance along them, so the cell's footprint
follows its transcripts rather than its segmentation mask alone.

Throughout, $d$ is the number of spatial dimensions, $2$ or $3$, and
$\Lambda_c = \Sigma_c^{-1}$ is the precision.

## The position term

The position factor in the [intensity](overview.md#the-generative-model) is a Gaussian
density, normaliser included:

$$
\begin{aligned}
e^{-D_c(x)} &= \mathcal{N}(x \mid m_c, \Sigma_c), \\[4pt]
D_c(x) &= \tfrac12 (x - m_c)^T \Lambda_c (x - m_c)
  + \tfrac12 \log|\Sigma_c| + \tfrac d2 \log 2\pi .
\end{aligned}
$$

Because it is a density, $A_c = \int e^{-D_c(x)}\,dx = 1$ whatever $m_c$ and $\Sigma_c$
are, with the inside-cell bonus off, so the expected-count term of the Poisson likelihood
does not depend on them. A Poisson
process with a Gaussian profile is the same as drawing a Poisson number of spots and then
placing each one independently at $x \sim \mathcal{N}(m_c, \Sigma_c)$. The cells are a
Gaussian mixture in space, with the mixture weights set by the expression model.

The only term of the [log-joint](overview.md#the-generative-model) that carries $m_c$ and
$\Sigma_c$ is therefore the log-intensity of the spots assigned to the cell:

$$
\log p(\dots) \supset \sum_s z_{s,c}\, \log \mathcal{N}(x_s \mid m_c, \Sigma_c) .
$$

The class indicator $\zeta_{c,k}$ multiplies this term in the log-joint, but
$\sum_k \zeta_{c,k} = 1$, so it sums out. The cell's class, scale factors and efficiency play
no part in where the cell sits.

## Prior

The conjugate prior for a Gaussian with unknown mean and covariance is the
Normal-Inverse-Wishart (NIW). The covariance gets an inverse Wishart prior and the mean,
conditional on the covariance, a Gaussian whose spread is that covariance shrunk by a
pseudo-count:

$$
p(m_c, \Sigma_c)
= \mathcal{N}\!\big(m_c \mid m_{0,c},\, \Sigma_c / \kappa_0\big)\;
  \mathcal{IW}\!\big(\Sigma_c \mid V_0,\, \nu_0\big),
$$

with

$$
\mathcal{IW}(\Sigma \mid V, \nu)
= \frac{|V|^{\nu/2}}{2^{\nu d/2}\,\Gamma_d(\nu/2)}\;
  |\Sigma|^{-(\nu + d + 1)/2}\,
  \exp\!\big(-\tfrac12 \operatorname{tr}(V \Sigma^{-1})\big),
$$

whose mean is $\mathbb{E}[\Sigma] = V / (\nu - d - 1)$ for $\nu > d + 1$.

The hyperparameters map onto the segmentation and two configuration settings:

- $m_{0,c}$ is the **segmentation centroid** of cell $c$;
- $\kappa_0$ is the **`cell_centroid_prior`** setting, the pseudo-count behind the centroid,
  $1000$ by default;
- $\nu_0$ is the **`cell_cov_prior`** setting, the degrees of freedom behind the covariance,
  $1000$ by default;
- $V_0 = (\nu_0 - d - 1)\, r^2 I$, with $r$ the mean cell radius, shared by every cell.

The scale matrix $V_0$ is chosen so that the prior mean of the covariance is exactly the
fixed covariance of the original model:

$$
\mathbb{E}[\Sigma_c] = \frac{V_0}{\nu_0 - d - 1} = r^2 I \quad\text{for any } \nu_0 > d + 1 .
$$

This is the parameterisation Zhang (2021) recommends. Setting $V_0 = \nu_0\, r^2 I$ instead,
which is the obvious first guess, puts the prior mean at $\tfrac{\nu_0}{\nu_0 - d - 1} r^2 I$
and makes the prior shape depend on $\nu_0$, so that changing the strength of the prior
would also change what it believes. With the choice above, $\nu_0$ only sets how strongly
$r^2 I$ is held, in the same way `rRho` holds `MisreadDensity` on the
[misread density](misread-density.md) page.

## Variational update

The [coordinate-ascent update](overview.md#the-variational-approximation) keeps the terms of
the log-joint that involve $m_c, \Sigma_c$ and takes the expectation over the rest. Only
the spot indicators $z_{s,c}$ carry any other factor, and
$\mathbb{E}[z_{s,c}] = \bar z_{s,c}$ is the spot's assignment probability from
[spot assignment](spot-assignment.md):

$$
\log q(m_c, \Sigma_c)
= \sum_s \bar z_{s,c}\, \log \mathcal{N}(x_s \mid m_c, \Sigma_c)
  + \log p(m_c, \Sigma_c) + \text{const} .
$$

Define the soft count, the soft mean and the soft scatter of the cell's spots:

$$
\bar n_c = \sum_s \bar z_{s,c},
\qquad
\bar x_c = \frac{1}{\bar n_c} \sum_s \bar z_{s,c}\, x_s,
\qquad
S_c = \sum_s \bar z_{s,c}\, (x_s - \bar x_c)(x_s - \bar x_c)^T .
$$

$\bar n_c$ is the cell's total expected count, summed over genes. Expanding the Gaussian and
splitting the quadratic about $\bar x_c$,

$$
\begin{aligned}
&\sum_s \bar z_{s,c}\, (x_s - m_c)^T \Lambda_c (x_s - m_c) \\[4pt]
&\qquad = \operatorname{tr}(\Lambda_c S_c)
  + \bar n_c\, (\bar x_c - m_c)^T \Lambda_c (\bar x_c - m_c),
\end{aligned}
$$

the data term becomes

$$
\begin{aligned}
&\sum_s \bar z_{s,c}\, \log \mathcal{N}(x_s \mid m_c, \Sigma_c) \\[4pt]
&\qquad = -\frac{\bar n_c}{2} \log|\Sigma_c|
  - \frac12 \operatorname{tr}(\Lambda_c S_c)
  - \frac{\bar n_c}{2} (\bar x_c - m_c)^T \Lambda_c (\bar x_c - m_c) + \text{const} .
\end{aligned}
$$

Adding the log prior,

$$
\begin{aligned}
\log p(m_c, \Sigma_c)
&= -\frac12 \log|\Sigma_c|
  - \frac{\kappa_0}{2} (m_c - m_{0,c})^T \Lambda_c (m_c - m_{0,c}) \\
&\quad - \frac{\nu_0 + d + 1}{2} \log|\Sigma_c|
  - \frac12 \operatorname{tr}(V_0 \Lambda_c) + \text{const},
\end{aligned}
$$

the two quadratics in $m_c$ combine into one. Completing the square with weights $\kappa_0$
on $m_{0,c}$ and $\bar n_c$ on $\bar x_c$ gives a Gaussian in $m_c$ with pseudo-count
$\kappa_0 + \bar n_c$, and leaves behind a cross term in $(\bar x_c - m_{0,c})$ that joins the
scale matrix. Collecting the powers of $|\Sigma_c|$ and the traces, the result is again a
Normal-Inverse-Wishart:

$$
\boxed{\;
q(m_c, \Sigma_c)
= \mathcal{N}\!\big(m_c \mid \hat m_c,\, \Sigma_c / \hat\kappa_c\big)\;
  \mathcal{IW}\!\big(\Sigma_c \mid \hat V_c,\, \hat\nu_c\big)
\;}
$$

with

$$
\begin{aligned}
\hat\kappa_c &= \kappa_0 + \bar n_c, \\[2pt]
\hat m_c &= \frac{\kappa_0\, m_{0,c} + \bar n_c\, \bar x_c}{\kappa_0 + \bar n_c}, \\[6pt]
\hat\nu_c &= \nu_0 + \bar n_c, \\[2pt]
\hat V_c &= V_0 + S_c
  + \frac{\kappa_0\, \bar n_c}{\kappa_0 + \bar n_c}\,
    (\bar x_c - m_{0,c})(\bar x_c - m_{0,c})^T .
\end{aligned}
$$

These are the conjugate updates of Bishop (2006, section 2.3.6) with the hard count $n$
replaced by the soft count $\bar n_c$ and the sample mean and scatter by their
$\bar z$-weighted versions, which is the form Bishop's section 10.2.1 uses for a mixture.
Note that
$S_c$ is the scatter about the soft mean $\bar x_c$, not about $\hat m_c$ or $m_{0,c}$; the
last term of $\hat V_c$ is what accounts for the centre having moved.

## Reading the result

**The centroid.** $\hat m_c$ is a weighted average of the segmentation centroid and the
soft mean of the cell's spots, with weights $\kappa_0$ and $\bar n_c$. A cell holding many
spots follows them; a cell holding few stays where the segmentation put it. $\kappa_0$ is
the number of spots the segmentation centroid counts for.

**The covariance.** The posterior mean is

$$
\begin{aligned}
\mathbb{E}[\Sigma_c] &= \frac{\hat V_c}{\hat\nu_c - d - 1} \\[6pt]
&= \frac{\bar n_c}{\bar n_c + \nu_0 - d - 1}\; \frac{\tilde S_c}{\bar n_c}
  + \frac{\nu_0 - d - 1}{\bar n_c + \nu_0 - d - 1}\; r^2 I ,
\end{aligned}
$$

where $\tilde S_c$ is $S_c$ with the centre-shift term folded in. This is equation (4) of
Zhang (2021): a weighted average of the spot scatter per spot, $\tilde S_c / \bar n_c$,
and the prior $r^2 I$, with weights $\bar n_c$ and $\nu_0 - d - 1$. So $\nu_0 - d - 1$, not $\nu_0$,
is the number of spots the prior shape counts for. With the default
$\nu_0 = 1000$ and $d = 3$ that is $996$; with $\nu_0 = 10$ it would be $6$.

$\hat V_c$ itself is a scatter, a sum of squared deviations, so its three parts add the way
counts add: $V_0$ is $\nu_0 - d - 1$ imaginary spots each worth $r^2$, and $S_c$ is
$\bar n_c$ real spots each worth their empirical covariance. The weighted average appears
on dividing by the total count.

**The centre-shift term.** $S_c$ measures the spread of the real spots about their own
mean $\bar x_c$, and the $\kappa_0$ imaginary spots behind the prior sit on the
segmentation centroid $m_{0,c}$ with no spread at all. But the cell's centre is $\hat m_c$,
in between. Measured from $\hat m_c$, the real spots spread by
$S_c + \bar n_c (\bar x_c - \hat m_c)(\bar x_c - \hat m_c)^T$ and the imaginary ones by
$\kappa_0 (m_{0,c} - \hat m_c)(m_{0,c} - \hat m_c)^T$. Substituting $\hat m_c$, the two
extra pieces sum to exactly the last term of $\hat V_c$,

$$
\begin{aligned}
&\bar n_c (\bar x_c - \hat m_c)(\bar x_c - \hat m_c)^T
+ \kappa_0 (m_{0,c} - \hat m_c)(m_{0,c} - \hat m_c)^T \\[4pt]
&\qquad = \frac{\kappa_0\, \bar n_c}{\kappa_0 + \bar n_c}\,
  (\bar x_c - m_{0,c})(\bar x_c - m_{0,c})^T .
\end{aligned}
$$

It is the parallel axis theorem: when two clouds are measured from a shared centre, each
picks up the squared distance from its own centre to the shared one. Two clouds that
disagree about where the centre is make a wider cloud together, and the cell comes out
stretched along the line from the segmentation centroid to the spots. The term vanishes
when the spots sit on the segmentation centroid, and when $\kappa_0 = 0$, where the
segmentation has no say and $\hat V_c = V_0 + S_c$.

**Empty cells.** A cell with $\bar n_c = 0$ has $\hat\kappa_c = \kappa_0$, $\hat m_c = m_{0,c}$,
$\hat\nu_c = \nu_0$ and $\hat V_c = V_0$: its posterior is its prior. No special case is
needed.

**The background** ($c = 0$) has no position and no shape, so there is no
$q(m_0, \Sigma_0)$. Its spots are scored by the [misread density](misread-density.md).

## The expected log density

The [spot assignment](spot-assignment.md) needs the expectation of the log density of a
spot's position under $q(m_c, \Sigma_c)$. Write it as $-\overline{D}_c(x_s)$ to mirror the
fixed $-D_c(x_s)$ of the original model:

$$
\begin{aligned}
-\overline{D}_c(x_s)
&= \mathbb{E}_{q}\big[\log \mathcal{N}(x_s \mid m_c, \Sigma_c)\big] \\[4pt]
&= -\frac d2 \log 2\pi
  + \frac12\, \mathbb{E}\big[\log|\Lambda_c|\big]
  - \frac12\, \mathbb{E}\big[(x_s - m_c)^T \Lambda_c (x_s - m_c)\big] .
\end{aligned}
$$

Under the NIW posterior, $\Lambda_c$ is Wishart with scale $\hat V_c^{-1}$ and $\hat\nu_c$
degrees of freedom, and $m_c \mid \Lambda_c$ is Gaussian about $\hat m_c$ with precision
$\hat\kappa_c \Lambda_c$. Both expectations are in closed form:

$$
\mathbb{E}\big[\log|\Lambda_c|\big]
= \sum_{i=1}^{d} \psi\!\left(\frac{\hat\nu_c + 1 - i}{2}\right)
  + d \log 2 - \log|\hat V_c| ,
$$

$$
\mathbb{E}\big[(x_s - m_c)^T \Lambda_c (x_s - m_c)\big]
= \frac{d}{\hat\kappa_c}
  + \hat\nu_c\, (x_s - \hat m_c)^T \hat V_c^{-1} (x_s - \hat m_c) ,
$$

with $\psi$ the digamma function (Bishop, 2006, equations 10.64 and 10.65, written there
for the Wishart on the precision). The second follows from averaging first over
$m_c$ given $\Lambda_c$, which adds $\operatorname{tr}(\Lambda_c \Sigma_c/\hat\kappa_c)
= d/\hat\kappa_c$, and then over $\Lambda_c$, using $\mathbb{E}[\Lambda_c] = \hat\nu_c \hat V_c^{-1}$.

So the position term is a quadratic in the spot's distance from $\hat m_c$ with precision
$\mathbb{E}[\Lambda_c] = \hat\nu_c \hat V_c^{-1}$, plus a constant for the cell. Two things
to note about it:

- It is **not** the log density at the posterior mean. Plugging $\hat m_c$ and
  $\mathbb{E}[\Sigma_c]$ into a Gaussian log density would use precision
  $\mathbb{E}[\Sigma_c]^{-1} = (\hat\nu_c - d - 1) \hat V_c^{-1}$, which is smaller than
  $\mathbb{E}[\Lambda_c]$ by the factor $(\hat\nu_c - d - 1)/\hat\nu_c$. As Zhang (2021)
  points out, the mean of the precision is not the inverse of the mean of the covariance.
  The coordinate-ascent update wants the former. This is the same distinction as
  $\overline{\log\rho}$ against $\log\bar\rho$ for the background.
- The term $d / (2\hat\kappa_c)$ is a flat cost for not knowing where the centre is. It is
  the same for every spot of a cell but differs between cells, so it tilts a spot towards
  the cell with more spots, and towards the background, which has no such term. It is at
  most $d / (2\kappa_0)$, the value for a cell with no spots: $0.0015$ in the log score at
  the default $\kappa_0 = 1000$ with $d = 3$. At $\kappa_0 = 10$ it is $0.15$ for a cell with
  no spots and $0.014$ for a cell with $100$.

## Contribution to the ELBO

The new factor adds $\mathbb{E}_q[\log p(m_c, \Sigma_c)] - \mathbb{E}_q[\log q(m_c, \Sigma_c)]
= -\mathrm{KL}\big(q(m_c, \Sigma_c)\,\|\,p(m_c, \Sigma_c)\big)$ per cell. The KL between two
NIW distributions splits into the conditional Gaussian part, averaged over $\Sigma_c$, and
the inverse Wishart part:

$$
\begin{aligned}
\mathrm{KL}_{\mathcal{N}}
&= \frac d2 \log\frac{\hat\kappa_c}{\kappa_0} - \frac d2 + \frac{d\,\kappa_0}{2\hat\kappa_c} \\
&\quad + \frac{\kappa_0\, \hat\nu_c}{2}\,
   (\hat m_c - m_{0,c})^T \hat V_c^{-1} (\hat m_c - m_{0,c}), \\[8pt]
\mathrm{KL}_{\mathcal{IW}}
&= \frac{\hat\nu_c}{2} \log|\hat V_c| - \frac{\nu_0}{2} \log|V_0|
   - \frac{(\hat\nu_c - \nu_0)\, d}{2} \log 2 \\
&\quad - \log\Gamma_d\!\Big(\frac{\hat\nu_c}{2}\Big)
   + \log\Gamma_d\!\Big(\frac{\nu_0}{2}\Big)
   + \frac{\hat\nu_c - \nu_0}{2}\, \mathbb{E}\big[\log|\Lambda_c|\big] \\
&\quad + \frac{\hat\nu_c}{2} \Big( \operatorname{tr}\big(V_0 \hat V_c^{-1}\big) - d \Big) .
\end{aligned}
$$

The position term inside the ELBO's data part is $-\overline{D}_c(x_s)$, the same expected log
density the assignment uses, so the ELBO and the updates agree.

## Simplifications

- **The update runs after the spot assignment.** The first spot assignment of a run is
  therefore made with the prior geometry, and every geometry update sees assignment
  probabilities that came out of the score. Run before it, the first update would see
  the uniform starting weights, under which the soft mean of a cell's spots is the
  centre of its whole neighbourhood, and the centres would jump there.
- **Nearest neighbours are not rebuilt.** Each spot competes only among its nearest cells,
  and that list is built once from the segmentation centroids. Once $\hat m_c$ moves, the
  list is slightly stale. It is a truncation of the sum over $c$ and is kept fixed for
  simplicity.
- **The inside-cell bonus keeps its fixed area factor.** With the bonus on, $A_c$ depends
  on the cell's shape and on the segmentation mask. The model treats it as the constant
  derived for a spherical cell of radius $r$, and ignores that the Gaussian may now be
  elongated or correlated across axes. With the bonus off, $A_c = 1$ exactly and nothing is
  ignored.

## References

Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*. Springer. Sections
2.3.6, 10.2 and appendix B.

Zhang, Z. (2021). A note on Wishart and inverse Wishart priors for covariance matrix.
*Journal of Behavioral Data Science*, 1(2), 119-126.
