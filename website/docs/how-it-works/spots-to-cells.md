---
id: spots-to-cells
title: Assigning spots to cells
sidebar_label: 4. Spots to cells
sidebar_position: 5
---

# Block 4: Assigning spots to cells

The final block of the loop assigns each RNA spot to the cell most likely to have
produced it, or to the background. It also closes the loop: once the spots are
reassigned, the gene counts per cell change, and the next iteration begins.

## The idea in one line

For each spot, evaluate the few nearest cells together with the background option, and
assign the spot to whichever explains it best.

## What determines the assignment

A spot is scored against each candidate cell using several terms, which combine
additively:

- **Distance.** Nearer cells are favoured. A spot located over a cell's nucleus is far
  more likely to belong to it than one several cell diameters away.

- **Compatibility with the gene.** If the cell's probable type expresses the spot's gene
  strongly, the spot is consistent with that cell; if the type rarely produces the gene,
  the spot is a poor match even when close. This is where the cell-type estimates from
  [block 3](cell-to-celltype.md) re-enter.

- **Cell-specific scaling.** The per-cell and per-gene-per-cell factors from
  [block 2](warping-the-reference.md) (theta and gamma) and the per-gene efficiency (eta)
  enter here, calibrating the expected fit so that the comparison is fair.

- **The inside-cell bonus.** A spot whose pixel lies within the cell's segmented boundary
  receives an additional term favouring that cell.

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