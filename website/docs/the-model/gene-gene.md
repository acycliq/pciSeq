
# Derivation: gene-gene dependence

The model so far keeps a structural feature of the original construction: gene counts
within a cell are conditionally independent given the class, so the
[cell-class update](cell-class.md) is a product $\prod_g \mathrm{NB}(\cdot)$ over genes.
That ignores co-regulation: if genes $A$ and $B$ are co-expressed in a class, a high count
of $A$ should raise the expectation of $B$ in the same cell. This extension relaxes the
independence assumption while preserving the closed-form conjugacy that drives cell typing.

## The latent shift and its covariance

We introduce a per-cell, gene-indexed latent vector $\mathbf{b}_c \in \mathbb{R}^G$ that
captures a correlated shift in expression across all genes. Conditional on the class, it
follows a centred multivariate normal whose covariance encodes the class-specific
co-expression structure:

$$
\mathbf{b}_c \mid k \sim \mathcal{N}(\mathbf{0}, \Sigma_k),
\qquad
\Sigma_k \sim \mathcal{IW}(\nu_0, \Psi_0) .
$$

The covariance $\Sigma_k$ itself is a random variable with a conjugate Inverse-Wishart
prior. The shift enters the intensity as a multiplicative, positive log-scale correction:

$$
\lambda_{g,c}(x) = \theta_c\, \mu_{g,k(c)}\, e^{b_{g,c}}\, e^{-D_c(x)}\, \gamma_{g,c}\, \eta_g .
$$

Because $\mathbf{b}_c$ is drawn from a covariance shared across all cells in a class, a
positive shift on gene $A$ correlates with a positive shift on any gene $B$ that
co-expresses with $A$. As with [$\theta_c$](scale-theta.md), treating both
$\mathbf{b}_c$ and $\gamma_{g,c}$ as full random variables would give an intractable
Poisson-Gamma-Lognormal mixture, so $\mathbf{b}_c$ is restricted to a class-conditional
point estimate $q(\mathbf{b}_c \mid k) = \delta(\mathbf{b}_c - \hat{\mathbf{b}}_{c\mid k})$
and obtained by MAP.

## Posterior of the class covariance $\Sigma_k$

Keeping only the terms of the log-joint that depend on $\Sigma_k$ (the multivariate-normal
prior on $\mathbf{b}_c$ and its own Inverse-Wishart prior), and using that the expectation
under the Dirac is exact, $\mathbb{E}[\mathbf{b}_c\mathbf{b}_c^{\mathrm{T}} \mid k]
= \hat{\mathbf{b}}_{c\mid k}\hat{\mathbf{b}}_{c\mid k}^{\mathrm{T}}$, the update is another
Inverse-Wishart:

$$
\boxed{\;
q(\Sigma_k) = \mathcal{IW}(\nu_k, \Psi_k)
\;}
$$

with parameters

$$
\nu_k = \nu_0 + \sum_c \bar\zeta_{c,k},
\qquad
\Psi_k = \Psi_0 + \sum_c \bar\zeta_{c,k}\, \hat{\mathbf{b}}_{c\mid k}\, \hat{\mathbf{b}}_{c\mid k}^{\mathrm{T}} .
$$

The degrees of freedom $\nu_k$ grow with the soft number of cells in class $k$, and the
scale matrix $\Psi_k$ accumulates a class-weighted scatter of the latent shifts, anchored
at the prior $\Psi_0$ - the same observed-over-expected flavour as the rest of the model.
The downstream MAP step needs the expected precision, which the Wishart mean gives in
closed form:

$$
\mathbb{E}_q[\Sigma_k^{-1}] = \nu_k\, \Psi_k^{-1} .
$$

## MAP estimate of the latent shift $\mathbf{b}_c$

Conditioning on the class and taking the expectation over the other factors, the objective
for $\mathbf{b}_c$ (with $\Lambda_{g,c\mid k} = \mu_{g,k}\, A_c\, \bar\gamma_{g,c}\, \bar\eta_g\, \bar\theta_c$) is

$$
\mathbb{E}[\mathcal{L}(\mathbf{b}_c \mid k)]
= - \sum_g \Lambda_{g,c\mid k}\, e^{b_{g,c}}
  + \sum_g \bar N_{c,g}\, b_{g,c}
  - \tfrac{1}{2}\, \mathbf{b}_c^{\mathrm{T}} (\nu_k \Psi_k^{-1})\, \mathbf{b}_c
  + \text{const}.
$$

Its gradient and Hessian are

$$
\nabla = \bar{\mathbf{N}}_c - \mathbf{\Lambda}_{c\mid k} \odot e^{\mathbf{b}_c} - (\nu_k \Psi_k^{-1})\,\mathbf{b}_c,
\qquad
\mathbf{H} = - \mathrm{diag}\big(\mathbf{\Lambda}_{c\mid k} \odot e^{\mathbf{b}_c}\big) - (\nu_k \Psi_k^{-1}),
$$

where $\odot$ is the element-wise product. The Hessian is the negation of a sum of two
positive-definite matrices, hence strictly negative definite, so the objective is strictly
concave with a unique maximum. We reach it by Newton-Raphson:

$$
\boxed{\;
\mathbf{b}_c^{(t+1)} = \mathbf{b}_c^{(t)} - \mathbf{H}^{-1}\, \nabla
\;}
$$

A single Newton step per outer variational iteration suffices: the variational loop is
already iterative, so the inner step need not converge inside each pass (a standard
Generalised EM construction). At the optimum each component $\hat b_{g,c}$ pushes the
predicted count toward the observed count, while the quadratic prior pulls it back toward
zero with strength $\nu_k \Psi_k^{-1}$; the off-diagonal entries of the precision propagate
evidence between correlated genes - the dependence the independent product could not
encode.

## How the other posteriors change

The factor $e^{\hat b_{g,c}}$ enters wherever the per-cell predicted gene count appeared.
The rates of $\gamma_{g,c}$, $\theta_c$, and $\eta_g$ each pick it up, and the
[cell-class update](cell-class.md) acquires both the $e^{\hat b_{g,c}}$ factor inside the
Negative Binomial mean and a multivariate-normal density term that penalises classes whose
covariance makes the cell's shift improbable:

$$
\begin{aligned}
q\big(k(c)=k\big) \propto\;
& \pi_k \exp\Big(\beta \sum_{c'\in\mathcal{N}_c} \bar\zeta_{c',k}\Big)
  \cdot \exp\Big(\mathbb{E}_q\big[\log p(\hat{\mathbf{b}}_{c\mid k} \mid \Sigma_k)\big]\Big) \\
& \cdot \prod_g \mathrm{NB}\big(\bar N_{c,g};\, r_\gamma,\, \mu_{g,k}\, A_c\, \bar\eta_g\, \bar\theta_c\, e^{\hat b_{g,c}}\big) .
\end{aligned}
$$

The [spot-to-cell assignment](spot-assignment.md) likewise gains a $\hat b_{g_s,c}$ term
inside its per-spot exponential. The standard identities behind these updates, and why the
Dirac construction is sound, are catalogued in the
[self-consistency appendix](appendix-self-consistency.md). The relation of this
construction to the published multivariate Negative Binomial literature is discussed in the
source document's appendix.
