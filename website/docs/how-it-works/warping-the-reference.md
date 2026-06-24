
# Block 2: Warping the cell type definitions

This block is the most subtle of the four, and the hardest to verify.

The cell type definitions are a lookup table: for every known cell type, they list the
average expression of every gene, measured in a **separate** scRNA-seq experiment. The
problem is that the in situ experiment in front of you does not behave exactly like that
separate experiment. Genes are detected at different efficiencies, some cells capture
more transcripts than others, and the overall scale is different. If you compared your
cells against the raw definitions, the match would be off for reasons that have nothing
to do with biology.

So before comparing, pciSeq **warps** the definitions to fit this experiment. It rescales
the expected expression numbers until they are on the same scale as what you actually
observe.

This calibration is also part of how the model separates signal from noise. Once the
technical losses are absorbed by these factors, genuine reads align with the adjusted
expected expression, while reads that match no calibrated cell type are left to be
accounted for as background. Putting the definitions on the right scale is therefore not
only about comparability; it is also what lets real expression be told apart from
technical noise.

## Why this block is harder to judge

The outputs of the other blocks can be inspected directly. Spot-to-cell assignments can
be overlaid on the image and assessed for spatial plausibility, and a cell's assigned
type can be compared against the expression of established marker genes. The adjustments
made in this block admit no comparable check. The inefficiency factors are **nuisance
parameters**: quantities the model must estimate in order to reach the results of
interest (the cell types and spot assignments), but which are not themselves reported.
They are latent, never observed, and there is no ground truth against which to validate
them. They are identified only indirectly, through the improvement they
produce in the agreement between cells and their assigned types. Their influence is
evident in the final result, but the adjustments themselves are not, which makes this
block intrinsically harder to validate than the assignment steps.

## The warp is a stack of scaling factors

The warp is not one number. It is a **family of correction factors**, each one rescaling
the expected expression at a different level of detail. pciSeq calls them
*inefficiencies*, because they mostly describe how much signal is lost relative to the
single-cell data.

From the broadest to the most specific:

- **Inefficiency** - a single constant applied to the whole reference, encoding the fact
  that in situ sequencing detects only a fraction (for example, around a fifth) of what
  scRNA-seq reports. It rescales **every gene in every cell** by this one factor, and is
  fully systemic: it does not distinguish between individual genes or cells.

- **eta** ($\eta_g$) - one factor **per gene**, shared across all cells. Some genes are
  detected more efficiently than others; eta captures that. It is the same for every
  cell, but different for every gene.

- **theta** ($\theta_{c,k}$) - one factor **per cell, for each candidate cell type**.
  Some cells simply yield more transcripts than the definitions predict, others fewer;
  theta is a single whole-cell **scalar** that stretches or shrinks that cell's expected
  counts across all its genes. It is worked out separately for every type the cell might be,
  because what counts as "expected" depends on which type you are testing it against.

- **gamma** ($\gamma_{g,c,k}$) - one factor **per gene, per cell, per candidate type**.
  This is the most fine-grained and idiosyncratic correction: it adjusts a single gene
  in a single cell, and again it is computed separately for each type that cell might be.
  It accounts for the residual mismatch that none of the broader factors can explain.

The four divide into two groups. The two broad factors, **Inefficiency** and **eta**, are
the same no matter what type a cell turns out to be. The two fine ones, **theta** and
**gamma**, are
**conditional on the class**: they are recomputed for each candidate type, because the
expectation they correct against is itself class-specific. This is why
[block 3](cell-to-celltype.md) can use them while it scores a cell against every type at
once.

## A pyramid of granularity

These factors can be arranged as a stack, ordered by how much of the experiment each one
covers. The broad, systemic factor sits at the base, affecting everything at once.
Higher up, the corrections get narrower and more specific, up to gamma at the apex, which
applies to just one gene, in one cell, under one candidate type.

<figure class="diagram">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 920 460" role="img" aria-label="A pyramid of the four scaling factors">
  <defs></defs>
  <line x1="70" y1="396" x2="70" y2="66" stroke="currentColor" stroke-opacity="0.5" stroke-width="1.5" />
  <path d="M 70,54 L 76,68 L 70,64 L 64,68 Z" fill="currentColor" fill-opacity="0.5" />
  <text class="ip-axis" transform="rotate(-90 90,230)" x="90" y="230" text-anchor="middle">GRANULARITY INCREASES</text>
  <path class="ip-tier" fill="#34d399" d="M 300,60.0 L 342.9,141.0 L 257.1,141.0 Z" />
  <path class="ip-tier" fill="#10b981" d="M 252.9,149.0 L 347.1,149.0 L 387.9,226.0 L 212.1,226.0 Z" />
  <path class="ip-tier" fill="#059669" d="M 207.9,234.0 L 392.1,234.0 L 432.9,311.0 L 167.1,311.0 Z" />
  <path class="ip-tier" fill="#047857" d="M 162.9,319.0 L 437.1,319.0 L 477.9,396.0 L 122.1,396.0 Z" />
  <line class="ip-leader" x1="343" y1="100.5" x2="540" y2="100.5" />
  <line class="ip-leader" x1="368" y1="187.5" x2="540" y2="187.5" />
  <line class="ip-leader" x1="413" y1="272.5" x2="540" y2="272.5" />
  <line class="ip-leader" x1="458" y1="357.5" x2="540" y2="357.5" />
  <g>
    <text x="550" y="96" class="ip-glyph">&#947; <tspan class="ip-name">(gamma)</tspan></text>
    <text x="550" y="116" class="ip-desc">scales gene g's count in cell c, per class k</text>
    <text x="550" y="183" class="ip-glyph">&#952; <tspan class="ip-name">(theta)</tspan></text>
    <text x="550" y="203" class="ip-desc">scales cell c's total count, per class k</text>
    <text x="550" y="268" class="ip-glyph">&#951; <tspan class="ip-name">(eta)</tspan></text>
    <text x="550" y="288" class="ip-desc">scales gene g's count across all cells</text>
    <text x="550" y="353" class="ip-glyph" font-size="20">Inefficiency</text>
    <text x="550" y="373" class="ip-desc">scales the whole experiment at once</text>
  </g>
</svg>
<figcaption>The same idea, four levels of granularity. Wide and systemic at the bottom, narrow and idiosyncratic at the top.</figcaption>
</figure>

Why have all four instead of one? Because mismatch occurs at all of these levels at once.
There is a global scale difference between the two technologies (handled at the base),
on top of that a per-gene detection pattern (eta), on top of that per-cell variation
(theta), and on top of all that, irreducible gene-by-cell noise (gamma). Each factor
absorbs the mismatch at its own scale, and the rest is left to the others.

## Inefficiencies: the common statistic

Although each acts at a different level of detail, **every inefficiency is the same
statistic**: a ratio of observed over expected.

$$
\text{factor} \;=\;
\frac{\text{what was actually observed}}{\text{what the model expected}}
$$

- **gamma** compares the observed counts of *one gene in one cell* against what the
  definitions predict for that same gene and cell, *assuming a given type*.
- **theta** compares the observed total counts of *one cell* against the total the
  definitions predict for it, *assuming a given type*.
- **eta** compares the observed counts of *one gene across all cells* against the total
  predicted for that gene.

In each case the factor is the observed quantity divided by the expected one, so it is
**greater than 1 when more was observed than predicted** and **less than 1 when less was
observed**. The only thing that differs between the factors is the level of aggregation
before the ratio is formed: a single gene-cell pair under one type, a whole cell under
one type, or a whole gene across all cells.

The block therefore reduces to a single principle applied at different scales: **observed
over expected.**

## What feeds in and what comes out

- **Feeds in:** the current gene counts per cell, the current cell-type estimates, and
  the raw cell type definitions.
- **Comes out:** a warped version of the expected expression, rescaled at every level,
  ready for [block 3](cell-to-celltype.md) to score cells against types.