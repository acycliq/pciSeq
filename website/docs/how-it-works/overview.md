
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

<figure class="diagram">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 500" role="img" aria-label="The pciSeq variational loop">
  <defs>
    <path id="vlTxtPath1" d="M 400,90 A 160,160 0 0,1 560,250" />
    <path id="vlTxtPath2" d="M 560,250 A 160,160 0 0,1 400,410" />
    <path id="vlTxtPath3" d="M 400,410 A 160,160 0 0,1 240,250" />
    <path id="vlTxtPath4" d="M 240,250 A 160,160 0 0,1 400,90" />
  </defs>
  <text x="400" y="255" text-anchor="middle" class="vl-center-title">THE VARIATIONAL LOOP</text>
  <circle cx="400" cy="250" r="100" fill="none" stroke="currentColor" stroke-width="0.5" stroke-dasharray="2 6" opacity="0.2" />
  <g>
    <path class="vl-arc-band" d="M 416.6,60.7 A 190,190 0 0,1 585.1,207.3 L 598.8,204.1 L 565.0,250.0 L 522.8,221.7 L 536.4,218.5 A 140,140 0 0,0 412.2,110.5 Z" />
    <text class="vl-text-label"><textPath href="#vlTxtPath1" startOffset="50%" text-anchor="middle">Misread Density</textPath></text>
    <path class="vl-arc-band vl-arc-highlight" d="M 589.3,266.6 A 190,190 0 0,1 442.7,435.1 L 445.9,448.8 L 400.0,415.0 L 428.3,372.8 L 431.5,386.4 A 140,140 0 0,0 539.5,262.2 Z" />
    <text class="vl-text-label"><textPath href="#vlTxtPath2" startOffset="50%" text-anchor="middle">Warping Reference</textPath></text>
    <path class="vl-arc-band" d="M 383.4,439.3 A 190,190 0 0,1 214.9,292.7 L 201.2,295.9 L 235.0,250.0 L 277.2,278.3 L 263.6,281.5 A 140,140 0 0,0 387.8,389.5 Z" />
    <text class="vl-text-label"><textPath href="#vlTxtPath3" startOffset="50%" text-anchor="middle">Cell Typing</textPath></text>
    <path class="vl-arc-band" d="M 210.7,233.4 A 190,190 0 0,1 357.3,64.9 L 354.1,51.2 L 400.0,85.0 L 371.7,127.2 L 368.5,113.6 A 140,140 0 0,0 260.5,237.8 Z" />
    <text class="vl-text-label"><textPath href="#vlTxtPath4" startOffset="50%" text-anchor="middle">Spot Assignment</textPath></text>
  </g>
</svg>
<figcaption>The variational loop. Each block feeds the next, and the last block feeds
back into the first. The loop runs until the spot assignments stop changing.</figcaption>
</figure>

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