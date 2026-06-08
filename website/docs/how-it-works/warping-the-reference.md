---
id: warping-the-reference
title: Warping the cell type definitions
sidebar_label: 2. Warping the definitions
sidebar_position: 3
---

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

![A pyramid of the four scaling factors](../../static/img/inefficiency-pyramid.svg)

<div className="docs-figure"><figcaption>The same idea, four levels of granularity. Wide and systemic at the bottom, narrow and idiosyncratic at the top.</figcaption></div>

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