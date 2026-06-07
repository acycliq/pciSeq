---
id: cell-to-celltype
title: Assigning cells to cell types
sidebar_label: 3. Cells to types
sidebar_position: 4
---

# Block 3: Assigning cells to cell types

By this point in the loop, two things are ready: the gene counts inside each cell, and a
warped reference that lives on the right scale for this experiment. Now pciSeq can ask
the question everyone cares about: **what type is each cell?**

## The idea in one line

Score every cell against every known cell type by how well the cell's gene counts match
the type's (warped) expected expression, then turn those scores into probabilities.

## How the scoring works

For a given cell and a given candidate type, pciSeq asks: *if this cell really were this
type, how surprised would I be by the counts I actually see?* A type whose expected
expression lines up with the cell's counts gets a high score; a type that predicts genes
the cell does not have, or misses genes it does have, gets a low one.

Each gene contributes its own piece of evidence, and the pieces add up across all genes.
A handful of strong marker genes can dominate the call, which is exactly what you want:
the presence or absence of a few characteristic transcripts is often what separates two
closely related types.

## Turning scores into probabilities

The raw scores are then passed through a **softmax**, which is just a way of converting a
list of scores into probabilities that add up to one. The result is, for every cell, a
probability spread across the types: maybe 80% one type, 15% a close relative, 5%
something else. That spread is informative on its own. A confident cell concentrates
almost all its probability on a single type; an ambiguous cell splits it.

## Two extra ingredients

Beyond the raw gene-expression match, pciSeq folds in two more things:

- **A prior.** Some cell types are simply more common than others. The prior nudges the
  scores so that, all else being equal, a cell is more likely to be called a common type
  than a rare one.

- **A spatial nudge (the MRF).** Cells of the same type tend to cluster together in
  tissue. pciSeq adds a small bonus when a cell's neighbours share a type, so isolated,
  biologically implausible calls are discouraged. This is a gentle pull toward spatial
  coherence, not a hard rule, and it is capped so it can never overwhelm the actual
  expression evidence.

## The "Zero" type

There is always one special class, called **Zero**, that expects no expression at all.
It is the home for cells that are mostly empty: debris, poorly segmented fragments, or
cells whose markers were not in the gene panel. Giving these somewhere to go keeps them
from forcing themselves onto a real type they do not belong to.

## What feeds in and what comes out

- **Feeds in:** gene counts per cell, the warped reference, the prior, and the
  neighbourhood layout.
- **Comes out:** a probability over cell types for every cell, used by
  [block 4](spots-to-cells.md) to judge which spots a cell would plausibly emit.