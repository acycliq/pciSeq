---
id: spots-to-cells
title: Assigning spots to cells
sidebar_label: 4. Spots to cells
sidebar_position: 5
---

# Block 4: Assigning spots to cells

The last block of the loop decides, for every RNA spot, **which cell it most likely came
from** - or whether it is just background noise. This is also the block that closes the
loop: once spots are reassigned, the gene counts per cell change, and the whole thing
starts over.

## The idea in one line

For each spot, weigh up the few nearby cells (plus the "background noise" option) and
hand the spot to whichever explains it best.

## What makes a cell a good home for a spot

A spot is scored against each candidate cell using several pieces of evidence, which add
together:

- **Distance.** Closer cells are favoured. A spot sitting right on top of a cell's
  nucleus is far more likely to belong to it than one several cell-widths away.

- **Does the cell express this gene?** If the cell's likely type is one that expresses
  this gene strongly, the spot fits comfortably. If the cell's type would almost never
  produce this gene, the spot fits poorly even when it is close by. This is where the
  cell-type guesses from [block 3](cell-to-celltype.md) feed back in.

- **The cell's own quirks.** The per-cell and per-gene-per-cell scaling factors from
  [block 2](warping-the-reference.md) (theta and gamma) and the per-gene efficiency
  (eta) all enter here, adjusting the expected fit so the comparison is fair.

- **The inside-cell bonus.** A spot whose pixel falls inside the cell's segmented
  boundary gets an extra push toward that cell.

## The background option

Every spot also competes against the **background noise** option from
[block 1](misread-density.md). If no nearby cell explains the spot better than pure
noise would, the spot is labelled background. This is how genuine misreads get filtered
out: they fail to beat the noise floor for any cell.

## Probabilities, not hard handovers

As with cell types, the scores go through a **softmax** and become probabilities. A spot
can be 90% one cell and 10% its neighbour, rather than forced entirely onto one. When
the gene counts are tallied up for the next round, a spot contributes its probability to
each cell, so an ambiguous spot is shared rather than gambled.

## Closing the loop

Reassigning the spots changes how many of each gene sit inside each cell. Those updated
counts are exactly what [block 1](misread-density.md) and
[block 2](warping-the-reference.md) need to start the next round. The loop tightens with
each pass, and when the spot probabilities stop moving, pciSeq has converged and returns
its answer.

## What feeds in and what comes out

- **Feeds in:** spot locations, cell-type probabilities, the warped reference and its
  scaling factors, and the per-gene background rates.
- **Comes out:** a probability over nearby cells (and background) for every spot, which
  becomes the updated gene counts that drive the next turn of the loop.