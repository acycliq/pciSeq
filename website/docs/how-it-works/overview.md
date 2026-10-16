# How it works

pciSeq treats every unknown as a latent variable of one Bayesian model: the cell of
origin of each spot, the class of each cell, the detection efficiency of each gene, and
the per-cell and per-gene scale factors. Each has a prior. The target of inference is
their joint posterior given the spots and the cell type definitions.

The posterior has no closed form. It is approximated by variational inference: the
approximating distribution is restricted to a mean-field family, one factor per group of
latent variables, and the member of the family closest to the posterior in
Kullback-Leibler divergence is selected. This is equivalent to maximising a lower bound
on the model evidence.

The optimal factors depend on one another, so they are fitted by coordinate ascent. Each
factor is updated to its optimum given the current values of the others, and the sweep
repeats until the estimates stop changing. The steps below are these updates.

## The variational loop

<figure class="diagram">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 500" role="img" aria-label="The pciSeq variational loop: misread density, warping the definitions, cell typing, the class weights when cell_type_prior is weighted, spot assignment and cell geometry, in that order">
  <defs>
    <path id="vlTxt0" d="M 412.6,90.5 A 160,160 0 0,1 522.6,147.2" />
    <path id="vlTxt1" d="M 530.2,359.3 A 170,170 0 0,0 550.8,171.5" />
    <path id="vlTxt2" d="M 394.1,419.9 A 170,170 0 0,0 503.5,384.9" />
    <path id="vlTxt3" d="M 266.0,354.7 A 170,170 0 0,0 357.4,414.6" />
    <path id="vlTxt4" d="M 249.9,170.2 A 170,170 0 0,0 246.6,323.2" />
    <path id="vlTxt5" d="M 278.3,146.1 A 160,160 0 0,1 377.7,91.6" />
  </defs>
  <text x="400" y="255" text-anchor="middle" class="vl-center-title">THE VARIATIONAL LOOP</text>
  <circle cx="400" cy="250" r="100" fill="none" stroke="currentColor" stroke-width="0.5" stroke-dasharray="2 6" opacity="0.2" />
  <g>
    <path class="vl-arc-band" d="M 411.6,60.4 A 190,190 0 0,1 539.0,120.4 L 549.2,110.9 L 539.9,162.6 L 492.2,164.1 L 502.4,154.5 A 140,140 0 0,0 408.5,110.3 Z" />
    <text class="vl-text-label"><textPath href="#vlTxt0" startOffset="50%" text-anchor="middle">Misread density</textPath></text>
    <path class="vl-arc-band" d="M 567.0,159.3 A 190,190 0 0,1 551.7,364.3 L 562.9,372.8 L 510.4,372.6 L 500.6,325.8 L 511.8,334.3 A 140,140 0 0,0 523.0,183.2 Z" />
    <text class="vl-text-label"><textPath href="#vlTxt1" startOffset="50%" text-anchor="middle">Warping definitions</textPath></text>
    <path class="vl-arc-band" d="M 518.3,398.7 A 190,190 0 0,1 403.3,440.0 L 403.6,454.0 L 371.3,412.5 L 402.2,376.0 L 402.4,390.0 A 140,140 0 0,0 487.2,359.6 Z" />
    <text class="vl-text-label"><textPath href="#vlTxt2" startOffset="50%" text-anchor="middle">Cell typing</textPath></text>
    <path class="vl-arc-band vl-arc-optional" d="M 355.6,434.8 A 190,190 0 0,1 256.6,374.7 L 246.0,383.8 L 257.1,332.5 L 304.9,332.7 L 294.3,341.8 A 140,140 0 0,0 367.3,386.1 Z" />
    <text class="vl-text-label vl-text-optional"><textPath href="#vlTxt3" startOffset="50%" text-anchor="middle">Class weights</textPath></text>
    <path class="vl-arc-band" d="M 230.0,334.8 A 190,190 0 0,1 227.8,169.7 L 215.1,163.8 L 266.5,153.0 L 285.8,196.8 L 273.1,190.8 A 140,140 0 0,0 274.7,312.5 Z" />
    <text class="vl-text-label"><textPath href="#vlTxt4" startOffset="50%" text-anchor="middle">Spot assignment</textPath></text>
    <path class="vl-arc-band" d="M 253.4,129.1 A 190,190 0 0,1 363.7,63.5 L 361.1,49.7 L 400.0,85.0 L 376.0,126.3 L 373.3,112.6 A 140,140 0 0,0 292.0,160.9 Z" />
    <text class="vl-text-label"><textPath href="#vlTxt5" startOffset="50%" text-anchor="middle">Cell geometry</textPath></text>
  </g>
</svg>
<figcaption>The updates, in the order they run. The class weights, dashed, are updated only with <code>cell_type_prior='weighted'</code>. The last feeds the first.</figcaption>
</figure>

1. **[Misread density.](misread-density.md)** The rate of background spots per gene,
   estimated from the spots currently assigned to the background.

2. **[Warping the cell type definitions.](warping-the-reference.md)** Per-gene detection
   efficiencies and per-cell and per-cell-gene scale factors that rescale the reference
   expression to this experiment. These are latent, with no observed counterpart, and
   are identified only through the fit between cells and classes.

3. **[Cell to class.](cell-to-celltype.md)** Every cell is scored against every class
   with the warped definitions, and the scores are normalised to probabilities.

   With `cell_type_prior='weighted'` the [class weights](../the-model/class-prior.md),
   the prior share of each class, are then updated from these probabilities. With the
   default, `'uniform'`, they stay equal and this step is skipped.

4. **[Spots to cells.](spots-to-cells.md)** Every spot is assigned to one of its
   neighbouring cells or to the background, as a probability, given the cell class
   probabilities.

5. **[Cell geometry.](../the-model/gaussian.md)** The centre and the covariance of each
   cell's Gaussian are updated from the spots now assigned to it. This comes after the
   spot assignment so that it sees real spot weights, never the uniform ones of the
   first pass.

The new spot assignments change the gene counts per cell, which enter step 1 of the next
sweep, and the new geometry enters the next spot assignment. A single pass would leave each step conditioned on the initial values of the
others; iterating propagates each update to the rest.

## Convergence

After each sweep the change in the spot-to-cell probabilities is measured as the maximum
absolute change over all spots and candidate cells. The loop stops when this falls below
`CellCallTolerance`, or at `max_iter`. The maximum is a strict criterion: one spot still
moving between two cells keeps the loop running after every other spot has settled,
where a mean would have stopped earlier.
