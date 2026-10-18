
# Derivation: the cell-class assignment $q(\zeta)$

This is the posterior over the type of each cell. It is also where the **spatial prior**
enters: the Markov Random Field term that rewards neighbouring cells for sharing a class.

## From a Poisson-Gamma mixture to a Negative Binomial

By the [CAVI update](overview.md#the-variational-approximation), the structured factor for the
class and its cell-gene scale is the expected
[log-joint](overview.md#the-generative-model) over the remaining latents:

$$
\log q(\zeta,\gamma) = \mathbb{E}_{z,\eta,\theta}\big[\log p(x, g, z, \zeta, \gamma, \eta, \theta)\big] + \text{const}.
$$

Because $\theta_c$ is a [point estimate](scale-factors.md#theta), its expectation is exact and
$\theta_c$ passes through as a constant. Taking the expectation and keeping the terms in
$\zeta$ and $\gamma$ leaves a Poisson-Gamma mixture in $\gamma$:

$$
\begin{aligned}
\log q(\zeta, \gamma)
= \sum_{c,k} \zeta_{c,k} \Big[ \sum_g \big[
& \underbrace{- \mu_{g,k}\, A_c\, \bar\eta_g\, \bar\theta_{c\mid k}\, \gamma_{g,c}
   + \bar N_{c,g}\log(\mu_{g,k}\, \bar\theta_{c\mid k}\, \gamma_{g,c})}_{\text{Poisson}} \\
& + \underbrace{(r_\gamma - 1)\log\gamma_{g,c} - r_\gamma\,\gamma_{g,c}}_{\text{prior: }\gamma \sim \mathrm{Gamma}(r_\gamma, r_\gamma)}
  \big] \\
& + \underbrace{\log\pi_k}_{\text{baseline log-prior}}
  + \underbrace{\beta \sum_{c'\in\mathcal{N}_c} \bar\zeta_{c',k}}_{\text{spatial log-prior}}
  + \text{const} \Big] .
\end{aligned}
$$

Integrating out the per-gene-per-cell factor $\gamma_{g,c}$ against its Gamma prior
collapses each gene's Poisson term into a Negative Binomial. The cell-class posterior is
therefore

$$
\boxed{\;
q\big(k(c)=k\big)
\propto
\Big(\prod_g \mathrm{NB}\big(\bar N_{c,g};\, r_\gamma,\, \mu_{g,k}\, A_c\, \bar\eta_g\, \bar\theta_{c\mid k}\big)\Big)
\cdot
\pi_k \exp\Big(\beta \sum_{c'\in\mathcal{N}_c} \bar\zeta_{c',k}\Big)
\;}
$$

The two factors are:

- a **Negative Binomial likelihood** over genes, with effective mean
  $\mu_{g,k}\, A_c\, \bar\eta_g\, \bar\theta_{c\mid k}$ - how well the cell's gene counts match the
  warped expectation for class $k$;
- the **baseline prior** $\pi_k$ scaled by the **spatial term**
  $\beta \sum_{c'\in\mathcal{N}_c} \bar\zeta_{c',k}$.

Under `cell_type_prior` set to `weighted`, $\pi$ is itself a random variable and the
$\log \pi_k$ in the score is its expectation $\mathbb{E}[\log \pi_k]$ under the
[Dirichlet posterior](class-prior.md).

## The spatial prior (MRF)

The spatial term comes from a Markov Random Field prior on the class assignments. In the
extended log-joint it appears as

$$
\beta \sum_{c,k} \sum_{c'\in\mathcal{N}_c} \mathbf{1}(\zeta_{c,k} = \zeta_{c',k} = 1)
= \beta \sum_{c,k} \zeta_{c,k} \sum_{c'\in\mathcal{N}_c} \zeta_{c',k},
$$

which induces the prior

$$
p(\zeta) \propto \exp\Big(\beta \sum_{c,k} \sum_{c'\in\mathcal{N}_c} \zeta_{c,k}\,\zeta_{c',k}\Big) .
$$

Here $\mathcal{N}_c$ is the set of nearest neighbours of cell $c$ (there are
`nNeighbors` of them, 9 by default) and $\beta \geq 0$ controls the strength of the
coupling. A cell is more readily assigned to a class its neighbours share. The term is
additive in log space, so it competes with the expression evidence rather than replacing
it.

### Weighting the neighbours by distance

The sum above weights all neighbours equally. With a weight per neighbour, so that a
touching cell counts more than one several diameters away, the support that multiplies
$\beta$ is

$$
S_{c,k} = \sum_{c' \in \mathcal{N}_c} w_{c,c'}\, \bar\zeta_{c',k} .
$$

The weight is a gaussian on the distance, with a width set per cell:

$$
w_{c,c'} \;\propto\; \exp\!\Big(-\big(d_{c,c'} / \sigma_c\big)^2\Big),
\qquad
\sigma_c = \operatorname{median}_{c' \in \mathcal{N}_c} d_{c,c'} ,
$$

then rescaled so that $\sum_{c'} w_{c,c'} = |\mathcal{N}_c|$. That normalisation keeps the
support on the same scale as the class probabilities it multiplies, so $\beta$ means the
same thing everywhere, and it puts a ceiling on the whole term: since $\bar\zeta \le 1$,

$$
S_{c,k} \;\le\; |\mathcal{N}_c| \qquad\text{and so}\qquad
\beta\, S_{c,k} \;\le\; \beta \cdot \texttt{nNeighbors} .
$$

$\sigma_c$ is taken from the cell's own neighbour distances, so closeness is relative to
the local cell density and one $\beta$ applies to dense and sparse regions alike.

**Example.** Cell $A$ has two neighbours, $B$ at distance 30 and $C$ at 60. The
implementation uses `nNeighbors` of them; two are shown here:

$$
\sigma_A = \operatorname{median}(30, 60) = 45, \qquad
w_B = e^{-(30/45)^2} = 0.641, \qquad
w_C = e^{-(60/45)^2} = 0.169 ,
$$

After rescaling to sum to 2, $w_B = 1.583$ and $w_C = 0.417$, so $B$ contributes to $A$'s
score with 3.8 times the weight of $C$.

### Pooling sister classes

Some classes in the reference are close relatives. `037 DG Glut` and `038 DG-PIR Ex IMN` are
an example: they share most of their markers, `037` is common in the dentate gyrus and `038`
is rare. A `038` cell usually sits inside a patch of `037` cells, so the MRF works against it.

**Example.** Take a cell whose 9 neighbours are 8 confident `037` cells and 1 confident `038`
cell. To keep the numbers simple, say all 9 are equally close (weight 1 each) and
$\beta = 1$. The neighbour bonus is then

- `037`: 8 (eight neighbours back it)
- `038`: 1 (one neighbour backs it)

Now say the cell's own reads fit `038` better: its log-likelihood is $-100$ under `038` and
$-103$ under `037`. Adding the bonus:

| | reads | bonus | total |
|---|---|---|---|
| `037` | $-103$ | $8$ | $-95$ |
| `038` | $-100$ | $1$ | $-99$ |

`037` wins by 4, and after the softmax the cell is `037` with probability 0.98. The
neighbours overrode the reads, and they carry no information here: both classes occupy the
same place, so a surround of `037` cells says nothing about whether this cell is `037` or
`038`.

`mrf_pooled_classes` addresses this. With the two classes in a group,

```python
opts = {
    "mrf_pooled_classes": [["037 DG Glut", "038 DG-PIR Ex IMN"]],
}
```

and a neighbour of either class backs both. The bonus becomes $8 + 1 = 9$ for each:

| | reads | bonus | total |
|---|---|---|---|
| `037` | $-103$ | $9$ | $-94$ |
| `038` | $-100$ | $9$ | $-91$ |

The bonus is the same for both so it cancels out, and the reads decide: `038` with
probability 0.95. The neighbours still count against every class outside the group, so a
patch of DG cells still lowers the probability of, for example, an interneuron class.

A group can hold more than two classes. Pooling does more than stop the neighbours choosing
inside the group: the votes across the group add up, and the total competes against every
class outside it. Take a cell whose 9 neighbours are 4 cells of a
pyramidal class `P` and 5 interneurons, each from a different interneuron class:

| | bonus for `P` | bonus for each interneuron class |
|---|---|---|
| no pooling | $4$ | $1$ |
| all interneurons in one group | $4$ | $5$ |

Without pooling the neighbours favour `P`. With the interneurons pooled they favour
"interneuron", although no single interneuron class has more than one neighbour. For a
pair like `037` and `038`, which occupy the same place, this does not matter. A large group
can win cells that none of its classes would have won alone, so groups should be kept small.

**In the maths.** The MRF prior rewards two neighbours for having the same class. Pooling
changes that to "the same class, or two classes in the same group". The support that
multiplies $\beta$ becomes

$$
\tilde S_{c,j} = \begin{cases}
\sum_{k \in G} S_{c,k} & \text{if class } j \text{ is in a group } G, \\
S_{c,j} & \text{otherwise,}
\end{cases}
$$

which in the example is $\tilde S_{037} = \tilde S_{038} = 8 + 1$. The ceiling from the
previous section still holds: a neighbour's class probabilities add up to at most 1, over a
group as well as over one class, so the term still never exceeds $\beta \cdot \texttt{nNeighbors}$.

**A class can only be in one group.**
