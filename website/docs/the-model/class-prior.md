
# Derivation: class prior $q(\pi)$

$\pi_k$ is the prior probability that a cell belongs to class $k$, the same for every cell.
$\pi$ without a subscript is the list of them, one entry per class, summing to one.
With `cell_type_prior` set to `uniform`, the default, $\pi$ is the fixed vector given by
`cell_type_weights`. With `weighted`, $\pi$ is a random variable with a Dirichlet prior,
and its posterior is refitted every iteration from the current class probabilities, so a
class that is common in the tissue gets a larger prior than a rare one. This page derives
that posterior.

The Zero class is held out. Its weight $w$ is fixed by `cell_type_weights`, and the Dirichlet
runs over the $K$ real classes, which share the remaining $1 - w$. From here on
$\pi = (\pi_1, \dots, \pi_K)$ is the list for the real classes only.

## Prior

A cell is Zero with probability $w$ and of real class $k$ with probability
$(1 - w)\,\pi_k$, where the $\pi_k$ sum to one. For cell $c$
with class indicator $\zeta_{c,k}$,

$$
\log p(\zeta_c \mid \pi)
= \zeta_{c,\text{Zero}} \log w
+ \sum_{k \ne \text{Zero}} \zeta_{c,k} \big[\log(1 - w) + \log \pi_k\big] .
$$

The prior on $\pi$ is a Dirichlet with pseudo-counts $\delta_k$,

$$
\log p(\pi) = \sum_{k \ne \text{Zero}} (\delta_k - 1) \log \pi_k + \text{const} .
$$

$\delta_k$ is the effective prior number of cells of class $k$ (Bishop 2006, section 10.2, p. 475).
It is fixed at $1$ for every real class and is not a setting. One cell per class is a very
weak prior against the thousands of cells of a dataset, so the class probabilities are set
by the data.

## Variational update

The terms of the log-joint that contain $\pi$ are the two above, the first summed over
cells. Taking
the expectation over every other factor replaces $\zeta_{c,k}$ by its posterior mean
$\bar\zeta_{c,k}$ and leaves

$$
\log q(\pi)
= \sum_{k \ne \text{Zero}} \Big( \sum_c \bar\zeta_{c,k} + \delta_k - 1 \Big) \log \pi_k + \text{const} .
$$

A Dirichlet has one parameter per class, and its log density is the sum over classes of
that parameter minus one, times $\log \pi_k$. The expression above has exactly that form,
so the posterior is a Dirichlet, like the prior. Its parameter for class $k$, written
$\alpha_k$, is what stands in front of $\log \pi_k$, plus one:

$$
\boxed{\;
q(\pi) = \mathrm{Dir}(\pi \mid \alpha_1, \dots, \alpha_K), \qquad
\alpha_k = \sum_c \bar\zeta_{c,k} + \delta_k
\;}
\tag{1}
$$

The prior had $\delta_k$ for class $k$; the posterior has $\alpha_k$, which is $\delta_k$
plus $\sum_c \bar\zeta_{c,k}$, the number of cells of class $k$ under the current
assignment. $\alpha$ without a subscript is the list of all the $\alpha_k$.

## Reading the result

The mean of a Dirichlet with parameters $\alpha$ is $\alpha_k / \sum_j \alpha_j$. With the
$\alpha_k$ of equation (1), the posterior mean of $\pi_k$ is

$$
\mathbb{E}[\pi_k] = \frac{\alpha_k}{\sum_j \alpha_j}
= \frac{\sum_c \bar\zeta_{c,k} + \delta_k}
{\sum_j \big[\sum_c \bar\zeta_{c,j} + \delta_j\big]} .
\tag{2}
$$

$\bar\zeta_{c,k}$ is the probability that cell $c$ has class $k$, so the sum
$\sum_c \bar\zeta_{c,k}$ is the number of cells of class $k$. Equation (2) is then the
share of the cells that are of class $k$: the more often a class appears, the greater its
weight in the model.

The prior and the assignment feed each other. A class that gains cells gets a larger
prior, which makes the next cell more likely to be assigned to it. Under `uniform` there
is no such feedback.

Including the Zero class, the prior probabilities are $w$ for Zero and
$(1 - w)\,\mathbb{E}[\pi_k]$ for a real class.

## Contribution to the class score

In the [cell-class update](cell-class.md), the term $\log \pi_k$ of the score of a real
class is replaced by its expectation under $q(\pi)$,

$$
\log(1 - w) + \mathbb{E}[\log \pi_k], \qquad
\mathbb{E}[\log \pi_k] = \psi(\alpha_k) - \psi\Big(\sum_j \alpha_j\Big),
$$

with $\psi$ the digamma function. For a class with more than a few cells (about ten) it is
close to $\log\big(\alpha_k / \sum_j \alpha_j\big)$, the log of the class's share of the cells. This uses

$$
\psi(x) \approx \log x - \frac{1}{2x} ,
$$

so the difference between $\psi$ and $\log$ is about $1/(2x)$: $0.05$ at $x = 10$, and
smaller beyond. For the Zero class the term is $\log w$. Its
weight is fixed and is not part of the Dirichlet.

$\alpha$ is refitted from the new class probabilities after every class update.

## References

Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*. Springer. Section 10.2,
the same update for the mixing coefficients of a Gaussian mixture.
