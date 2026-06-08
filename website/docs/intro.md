---
id: intro
title: What is pciSeq?
sidebar_label: Introduction
sidebar_position: 1
slug: /
---

# What is pciSeq?

**pciSeq** (probabilistic cell typing by in situ sequencing) processes the output of an
in situ sequencing experiment and addresses two coupled questions simultaneously:

1. **Which cell does each RNA spot belong to?**
2. **What is the type of each cell?**

The two questions are interdependent. The type of a cell cannot be determined with
confidence without knowing which spots lie within it, and a spot cannot be assigned to a
cell without an estimate of that cell's type. pciSeq resolves this by estimating both
quantities jointly, refining them in alternation until they stabilise.

## Inputs

- **Spots.** A table of detected RNA reads, each with its gene identity and spatial
  location (`x`, `y`, and a `z`-plane for 3D data).
- **A segmentation.** A label image indicating which pixels belong to which cell,
  typically derived from a DAPI nuclear stain.
- **Cell type definitions.** A table of average expression per gene for each known cell
  type, obtained from a separate scRNA-seq experiment. These provide the reference
  profiles of the candidate cell types.

## Outputs

- A **cell type** for every cell, expressed as a probability distribution over the known
  types, which also conveys the confidence of the assignment.
- A **parent cell** for every spot, again as a probability: a spot may be shared between
  neighbouring cells or attributed to the background.

## Why probabilities rather than hard labels

Segmentation boundaries are imprecise, gene detection is imperfect, and a fraction of
reads are noise. Rather than committing to single answers, pciSeq represents every
assignment as a probability and lets the evidence accumulate across iterations. A spot
that lies clearly within one cell and matches its expression is assigned with high
confidence; an ambiguous spot near a boundary is divided between cells. The same applies
to cell types.

## Where to go next

The following pages describe how the algorithm works, one building block at a time. Begin
with the [overview](how-it-works/overview.md) for the structure of the loop, then read
the blocks in order.