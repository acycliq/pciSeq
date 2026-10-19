
# Derivation: cell geometry $q(m_c, \Sigma_c)$

Each cell has a Gaussian with centre $m_c$ and covariance $\Sigma_c$, and a spot is scored
against the cell by the density of its position under it. In the original model the centre
is the segmentation centroid and the covariance is $r^2 I$, with $r$ the cell radius, set by `cell_radius`, by default half the
mean equivalent radius $\sqrt{\text{area}/\pi}$ of the segmented cells.
Here $m_c$ and $\Sigma_c$ are random variables, one pair per cell, with a prior and a
variational posterior. The Gaussian can then become an ellipsoid, elongated along the spots
of the cell, and its centre can shift where the segmentation did not capture the cell well.

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
\log p(\dots) = \sum_s z_{s,c}\, \log \mathcal{N}(x_s \mid m_c, \Sigma_c) + \text{const} .
$$

The class indicator $\zeta_{c,k}$ multiplies this term in the log-joint, but
$\sum_k \zeta_{c,k} = 1$, so it sums out. The cell's class, scale factors and efficiency do
not enter the update of its centre and covariance.

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
- $\kappa_0$ is the **`cell_centroid_prior`** setting, the pseudo-count of the prior on the centroid,
  $1000$ by default;
- $\nu_0$ is the **`cell_cov_prior`** setting, the degrees of freedom of the prior on the covariance,
  $1000$ by default;
- $V_0 = (\nu_0 - d - 1)\, r^2 I$, with $r$ the cell radius, shared by every cell.

The scale matrix $V_0$ is chosen so that the prior mean of the covariance is exactly the
fixed covariance of the original model:

$$
\mathbb{E}[\Sigma_c] = \frac{V_0}{\nu_0 - d - 1} = r^2 I \quad\text{for any } \nu_0 > d + 1 .
$$

This is the parameterisation Zhang (2021) recommends. Setting $V_0 = \nu_0\, r^2 I$ instead
puts the prior mean at $\tfrac{\nu_0}{\nu_0 - d - 1} r^2 I$
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
\bar N_c = \sum_s \bar z_{s,c},
\qquad
\bar x_c = \frac{1}{\bar N_c} \sum_s \bar z_{s,c}\, x_s,
\qquad
S_c = \sum_s \bar z_{s,c}\, (x_s - \bar x_c)(x_s - \bar x_c)^T .
$$

$\bar N_c$ is the cell's total expected count, summed over genes. Expanding the Gaussian and
splitting the quadratic about $\bar x_c$,

$$
\begin{aligned}
&\sum_s \bar z_{s,c}\, (x_s - m_c)^T \Lambda_c (x_s - m_c) \\[4pt]
&\qquad = \operatorname{tr}(\Lambda_c S_c)
  + \bar N_c\, (\bar x_c - m_c)^T \Lambda_c (\bar x_c - m_c),
\end{aligned}
$$

the data term becomes

$$
\begin{aligned}
&\sum_s \bar z_{s,c}\, \log \mathcal{N}(x_s \mid m_c, \Sigma_c) \\[4pt]
&\qquad = -\frac{\bar N_c}{2} \log|\Sigma_c|
  - \frac12 \operatorname{tr}(\Lambda_c S_c)
  - \frac{\bar N_c}{2} (\bar x_c - m_c)^T \Lambda_c (\bar x_c - m_c) + \text{const} .
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
on $m_{0,c}$ and $\bar N_c$ on $\bar x_c$ gives a Gaussian in $m_c$ with pseudo-count
$\kappa_0 + \bar N_c$, and leaves behind a cross term in $(\bar x_c - m_{0,c})$ that joins the
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
\hat\kappa_c = \kappa_0 + \bar N_c
\tag{1}
$$

$$
\hat m_c = \frac{\kappa_0\, m_{0,c} + \bar N_c\, \bar x_c}{\kappa_0 + \bar N_c}
\tag{2}
$$

$$
\hat\nu_c = \nu_0 + \bar N_c
\tag{3}
$$

$$
\hat V_c = V_0 + S_c
  + \frac{\kappa_0\, \bar N_c}{\kappa_0 + \bar N_c}\,
    (\bar x_c - m_{0,c})(\bar x_c - m_{0,c})^T
\tag{4}
$$

These are the conjugate updates of Bishop (2006, section 2.3.6) with the hard count $n$
replaced by the soft count $\bar N_c$ and the sample mean and scatter by their
$\bar z$-weighted versions, which is the form Bishop's section 10.2.1 uses for a mixture.
$S_c$ is the scatter about the soft mean $\bar x_c$, not about $\hat m_c$ or $m_{0,c}$; the
last term of $\hat V_c$ is what accounts for the centre having moved.

## Reading the result

**The fitted centroid.** Equation (2) is a weighted average of the segmentation centroid and the
soft mean of the cell's spots, with weights $\kappa_0$ and $\bar N_c$: $\kappa_0$ is the
number of spots the segmentation centroid counts for. For $\kappa_0 \gg \bar N_c$,
$\hat m_c \to m_{0,c}$.

**The covariance.** Write the last two terms of equation (4) as one,

$$
\tilde S_c = S_c
  + \frac{\kappa_0\, \bar N_c}{\kappa_0 + \bar N_c}\,
    (\bar x_c - m_{0,c})(\bar x_c - m_{0,c})^T ,
$$

so that $\hat V_c = V_0 + \tilde S_c$. The mean of an inverse Wishart is
$V / (\nu - d - 1)$. With $V_0 = (\nu_0 - d - 1)\, r^2 I$ from the prior and
$\hat\nu_c = \nu_0 + \bar N_c$ from equation (3),

$$
\begin{aligned}
\mathbb{E}[\Sigma_c] &= \frac{\hat V_c}{\hat\nu_c - d - 1} \\[6pt]
&= \frac{(\nu_0 - d - 1)\, r^2 I + \tilde S_c}{\bar N_c + \nu_0 - d - 1} \\[6pt]
&= \frac{\bar N_c}{\bar N_c + \nu_0 - d - 1}\; \frac{\tilde S_c}{\bar N_c}
  + \frac{\nu_0 - d - 1}{\bar N_c + \nu_0 - d - 1}\; r^2 I .
\end{aligned}
\tag{5}
$$

The last line multiplies and divides $\tilde S_c$ by $\bar N_c$. This is equation (4) of
Zhang (2021): a weighted average of the spot scatter per spot, $\tilde S_c / \bar N_c$,
and the prior $r^2 I$, with weights $\bar N_c$ and $\nu_0 - d - 1$. The pseudo-count of
the prior shape is $\nu_0 - d - 1$, not $\nu_0$: `cell_cov_prior = 10` in 3D counts as $6$
spots, and the setting must be above $d + 1$.

For $\nu_0 \gg \bar N_c$, $\mathbb{E}[\Sigma_c] \to r^2 I$.

**The centre-shift term.** $S_c$ is the scatter of the spots about their own mean
$\bar x_c$, but the fitted centroid of the cell is $\hat m_c$. About $\hat m_c$ the scatter of the
spots is larger by $\bar N_c (\bar x_c - \hat m_c)(\bar x_c - \hat m_c)^T$, and the prior
centre $m_{0,c}$, with pseudo-count $\kappa_0$, adds
$\kappa_0 (m_{0,c} - \hat m_c)(m_{0,c} - \hat m_c)^T$. Substituting $\hat m_c$ from
equation (2), the two sum to the last term of equation (4):

$$
\begin{aligned}
&\bar N_c (\bar x_c - \hat m_c)(\bar x_c - \hat m_c)^T
+ \kappa_0 (m_{0,c} - \hat m_c)(m_{0,c} - \hat m_c)^T \\[4pt]
&\qquad = \frac{\kappa_0\, \bar N_c}{\kappa_0 + \bar N_c}\,
  (\bar x_c - m_{0,c})(\bar x_c - m_{0,c})^T .
\end{aligned}
$$

The term stretches the covariance along the line from the segmentation centroid to the mean
of the spots. It is zero when $\bar x_c = m_{0,c}$, and when $\kappa_0 = 0$, where
$\hat V_c = V_0 + S_c$.

**Empty cells.** A cell with $\bar N_c = 0$ has $\hat\kappa_c = \kappa_0$, $\hat m_c = m_{0,c}$,
$\hat\nu_c = \nu_0$ and $\hat V_c = V_0$: its posterior is its prior. No special case is
needed.

**The background** ($c = 0$) has no position and no shape, so it has no geometry
posterior. Its spots are scored by the [misread density](misread-density.md).

## Choosing the prior strengths

Both settings are $1000$ by default. Each is a pseudo-count: the prior on the centroid counts
as $1000$ spots placed at the segmentation centroid, and the prior on the shape as about
$1000$ spots spread as a sphere of radius $r$. Cells
usually hold far fewer spots than that, so at the defaults the prior tends to dominate: the
centre stays close to the segmentation centroid and the shape close to the sphere.

Run with the defaults first. If the results are not satisfactory, lower the settings to the
order of the number of spots in a typical cell, or below. The two are independent:
`cell_centroid_prior` controls the centre, `cell_cov_prior` the shape.

Lower values make the run take more iterations: the spots a cell holds depend on its shape, and its shape depends on the
spots it holds. How many depends on the data and on the
other settings. The change between iterations may not fall steadily: it can remain at a
similar value for many iterations before it falls again, so a larger `max_iter` may be
needed. A cell can also widen and be assigned spots that fixed shapes assign elsewhere; the spot
assignment, in particular the fraction of background spots, should therefore be compared
with a run at the defaults.

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
$\mathbb{E}[\Lambda_c] = \hat\nu_c \hat V_c^{-1}$, plus a constant for the cell. It has two
properties:

- It is **not** the log density at the posterior mean. Plugging $\hat m_c$ and
  $\mathbb{E}[\Sigma_c]$ into a Gaussian log density would use precision
  $\mathbb{E}[\Sigma_c]^{-1} = (\hat\nu_c - d - 1) \hat V_c^{-1}$, which is smaller than
  $\mathbb{E}[\Lambda_c]$ by the factor $(\hat\nu_c - d - 1)/\hat\nu_c$. As Zhang (2021)
  points out, the mean of the precision is not the inverse of the mean of the covariance.
  The coordinate-ascent update requires the former. This is the same distinction as
  $\overline{\log\rho}$ against $\log\bar\rho$ for the background.
- The term $d / (2\hat\kappa_c)$ is a flat cost for not knowing where the centre is. It is
  the same for every spot of a cell but differs between cells, so it favours
  the cell with more spots, and the background, which has no such term. It is at
  most $d / (2\kappa_0)$, the value for a cell with no spots: $0.0015$ in the log score at
  the default $\kappa_0 = 1000$ with $d = 3$. At $\kappa_0 = 10$ it is $0.15$ for a cell with
  no spots and $0.014$ for a cell with $100$.

## Simplifications

- **The update runs after the spot assignment.** The first spot assignment of a run is
  therefore made with the prior geometry, and every geometry update sees assignment
  probabilities that came out of the score. Run before it, the first update would see
  the uniform starting weights, under which the soft mean of a cell's spots is the
  centre of its whole neighbourhood, and the centres would jump there.
- **Nearest neighbours are not rebuilt.** Each spot is scored only against its nearest cells,
  and that list is built once from the segmentation centroids. Once $\hat m_c$ moves, the
  list no longer matches the fitted centroids exactly. It is a truncation of the sum over $c$ and is kept fixed for
  simplicity.
- **The inside-cell bonus keeps its fixed area factor.** With the bonus on, $A_c$ depends
  on the cell's shape and on the segmentation mask. The model treats it as the constant
  computed at the start from the cell's segmented area, as for a round cell, and ignores that the Gaussian may now be
  elongated or correlated across axes. With the bonus off, $A_c = 1$ exactly and nothing is
  ignored.

## References

Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*. Springer. Sections
2.3.6, 10.2 and appendix B.

Zhang, Z. (2021). A note on Wishart and inverse Wishart priors for covariance matrix.
*Journal of Behavioral Data Science*, 1(2), 119-126.
