
# Block 4: Assigning spots to cells

The final block of the loop assigns each RNA spot to the cell most likely to have
produced it, or to the background. It also closes the loop: once the spots are
reassigned, the gene counts per cell change, and the next iteration begins.

## The idea in one line

For each spot, evaluate the few nearest cells together with the background option, and
assign the spot to whichever explains it best.

## What determines the assignment

When a cell weighs up a spot, it asks two questions: **where are you?** and **what are
you?** The first is geometry; the second is about identity - whether a cell like this would
produce this gene. The score adds the two together.

<figure class="diagram">
<img src="/img/spot-assignment-blocks.svg" alt="The building blocks of the spot-to-cell score">
<figcaption>The score, block by block. One block asks <em>where</em> the spot is; four ask <em>what</em> it is. The score simply adds them up.</figcaption>
</figure>

### Where: the spatial fit

A spot near a cell's centre scores higher than one several cell-widths away. Each cell has a
Gaussian shape, and this term is simply how well the spot's position sits inside that shape.
It is a hard geometric measurement, and it says nothing about which gene the spot carries.

### What: four facets

The "what" splits into four, each a different way of asking *does this gene belong in this
cell?* All four are weighted by how **confident** we are about the cell's type (from
[block 3](cell-to-celltype.md)): the surer the cell is of what it is, the more decisively
each one speaks.

**Alignment - does the cell's *type* express this gene?** If the cell is probably a type
that makes the gene strongly, the spot fits; if its type rarely produces the gene, the spot
is a poor match even when it sits right on the cell. Picture a spot midway between two cells,
all else equal: it is drawn to the one whose likely type expresses the gene. That overlap -
between what the cell probably is and what the gene marks - is its *alignment*.

**Gravity - is this a big, active cell?** Some cells gather more transcripts than their type
predicts (the per-cell scaling from [block 2](warping-the-reference.md)). Read that as the
cell's **mass**: a heavier cell pulls harder, so with everything else equal a spot drifts
toward whichever cell is already capturing the most. A "rich-get-richer" pull that lets a
clearly active cell claim the ambiguous spots around it.

**Enrichment - does *this* cell already carry this gene?** Easy to confuse with alignment,
but they differ. Alignment is about the cell's **type**, and is the same for every cell of
that type. Enrichment is about **this individual cell's own data** - whether it carries more
of the gene than its type predicts. The clean test: two cells of the *same* type have
identical alignment, so alignment cannot choose between them; but the cell already loaded
with the gene has the higher enrichment, and it wins the spot. **Alignment tells different
types apart; enrichment tells same-type cells apart.**

**Misread correction - is the gene reliably detected?** Some genes are read out far more
efficiently than others. This term does not choose between cells - it is the same for all of
them - but it matters in the contest against the **background** below, attenuating a poorly
detected gene's signal so its spots are more readily called misreads.

An optional **inside-cell bonus** can be switched on to additionally favour a spot whose
pixel falls within a cell's segmented boundary; it is off by default.

## The background option

Each spot also competes against the background option from
[block 1](misread-density.md). If no nearby cell explains the spot better than the
background, the spot is attributed to the background. This is how genuine misreads are
filtered out: they fail to exceed the background level for any cell.

## Probabilities, not hard assignments

As with cell types, the scores are passed through a **softmax** and become probabilities.
A spot may be assigned 0.9 to one cell and 0.1 to a neighbour rather than entirely to
one. When the gene counts are accumulated for the next iteration, each spot contributes
its probability to each cell, so an ambiguous spot is shared rather than committed to a
single cell.

## Closing the loop

Reassigning the spots changes how many copies of each gene fall within each cell. These
updated counts are the inputs that [block 1](misread-density.md) and
[block 2](warping-the-reference.md) require for the next iteration. The estimates are
refined on each pass, and when the spot probabilities stop changing the algorithm has
converged and returns its result.

## What feeds in and what comes out

- **Feeds in:** spot locations, cell-type probabilities, the warped definitions and their
  scaling factors, and the per-gene background rates.
- **Comes out:** a probability distribution over the nearby cells (and the background) for
  every spot, which becomes the updated gene counts that drive the next iteration.