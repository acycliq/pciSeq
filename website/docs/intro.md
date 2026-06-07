---
id: intro
title: What is pciSeq?
sidebar_label: Introduction
sidebar_position: 1
slug: /
---

# What is pciSeq?

**pciSeq** (probabilistic cell typing by in situ sequencing) takes the raw output of
an in situ sequencing experiment and answers two questions at the same time:

1. **Which cell does each RNA spot belong to?**
2. **What cell type is each cell?**

These two questions are tangled together. You cannot confidently say what type a cell
is until you know which spots are inside it, and you cannot confidently assign a spot
to a cell until you have a sense of what that cell is. pciSeq deals with this by
**going back and forth** between the two answers, refining both a little at a time
until they stop changing. That back-and-forth is the heart of the method.

## What goes in

- **Spots.** A table of detected RNA reads: which gene each one is, and where it sits
  in space (`x`, `y`, and a `z`-plane for 3D data).
- **A segmentation.** A label image marking which pixels belong to which cell, usually
  from a DAPI nuclear stain.
- **A single-cell reference.** A table of average expression per gene for each known
  cell type, taken from a separate scRNA-seq experiment. This is the "dictionary" of
  what each cell type is supposed to look like.

## What comes out

- A **cell type** for every cell, given as a probability over the known types (so you
  also see how confident the call is).
- A **parent cell** for every spot, again as a probability (a spot can be shared
  between neighbours, or flagged as background noise).

## Why probabilities, not hard labels

Segmentation boundaries are fuzzy, genes are detected imperfectly, and some reads are
just noise. Instead of pretending these are certain, pciSeq keeps everything as a
probability and lets the evidence accumulate. A spot that clearly sits inside one cell
and matches its expression gets assigned with high confidence; an ambiguous spot near a
boundary is split. The same is true for cell types.

## Where to go next

The rest of these docs walk through **how the algorithm actually works**, one building
block at a time. Start with the [overview](how-it-works/overview.md) to see the whole
loop, then read each block in order.