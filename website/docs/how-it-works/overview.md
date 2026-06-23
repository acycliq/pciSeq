
# How it works: the big picture

pciSeq does not compute its result in a single pass. It is an iterative procedure: a
small set of steps is applied repeatedly, each step refining its own estimate from the
current estimates of the others, until the estimates stabilise. The converged state is
the output.

pciSeq is a **fully Bayesian** model. Every unknown in the problem - the cell of origin
of each spot, the type of each cell, the gene detection efficiencies, and the per-cell
and per-gene scale factors - is treated as a latent variable with a prior distribution.
The target of inference is their **joint posterior**: the distribution over all of these
latent variables together, conditioned on the observed spots and the cell type
definitions.

This posterior is analytically intractable and cannot be evaluated in closed form.
**Variational inference** (variational Bayes) addresses this by approximation. We
restrict attention to a tractable family of distributions, factorised into one factor
per group of latent variables (a mean-field approximation), and select the member of
that family closest to the true posterior, where closeness is measured by the
Kullback-Leibler divergence. Equivalently, this maximises a lower bound on the model
evidence.

The optimal factors are coupled and cannot be obtained simultaneously, so the
approximation is fitted by **coordinate ascent**: an iterative loop that updates one
factor at a time, each to its optimal form given the current estimates of all the
others. Each latent variable is therefore estimated **conditionally on the rest**. Every
sweep tightens the approximation, and the loop runs until the estimates converge. The
four blocks below are exactly these conditional updates.

## The four building blocks

![The pciSeq variational loop](/img/variational-loop.svg)

<div className="docs-figure">
<figcaption>The variational loop. Each block feeds the next, and the last block feeds
back into the first. The loop runs until the spot assignments stop changing.</figcaption>
</div>

1. **[Estimate the misread density per gene.](misread-density.md)**
    Estimate how much background noise each gene produces, so genuine signal can be
   separated from it.

2. **[Warp the cell type definitions.](warping-the-reference.md)**
   Rescale the cell type definitions so they match the scale and characteristics of
   *this* experiment. This is the most subtle block, and the hardest to verify, because
   it happens entirely behind the scenes.

3. **[Assign cells to cell types.](cell-to-celltype.md)**
   With the warped definitions in hand, score every cell against every known type and
   turn the scores into probabilities.

4. **[Assign spots to cells.](spots-to-cells.md)**
   With cell types in hand, decide which cell each spot most likely came from (or
   whether it is background noise).

Then the loop closes: new spot assignments change the gene counts per cell, which feeds
straight back into block 1, and the cycle repeats.

## Why a loop and not a pipeline

Running the steps only once, in sequence, would leave each one based on crude initial
estimates of the others: the misread density would rest on a provisional spot assignment,
the cell types on uncalibrated definitions, and so on. Iterating allows each correction
to propagate. An improved misread estimate refines the spot assignments, which refine the
per-cell gene counts, which refine the cell-type estimates, which in turn refine the spot
assignments. Successive passes continue until the estimates no longer change
appreciably.

## Convergence

After each iteration, pciSeq measures how much the **spot-to-cell probabilities** have
changed. When this change falls below a fixed tolerance, the estimates are taken to have
converged and the loop terminates. A maximum number of iterations is also imposed as a
safeguard.

## A note on block 2

Three of the four blocks produce quantities that can be inspected directly: background
rates, cell-type scores, and spot assignments. Block 2 is different: the warped
definitions it produces are fully latent, with no observed counterpart, and are
identified only through their effect on the agreement between cells and types. It is where
the [family of inefficiency factors](warping-the-reference.md) is estimated.