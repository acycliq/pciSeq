---
id: misread-density
title: Estimating the misread density
sidebar_label: 1. Misread density
sidebar_position: 2
---

# Block 1: Estimating the misread density

Not every detected spot is real. Some reads are technical artefacts, optical crosstalk,
or genes decoded by mistake. pciSeq needs a way to say "this spot looks more like noise
than signal" and quietly set it aside. The misread density is how it does that.

## The idea in one line

For each gene, estimate **how many of its spots are just background noise spread evenly
across the tissue**, and use that as the bar a real spot has to clear.

## How it works

Imagine sprinkling noise spots uniformly over the whole tissue section. A gene with a
high background rate will scatter spots everywhere, even in empty space far from any
cell. A gene with a low background rate will almost never do that.

When pciSeq later asks "did this spot come from a nearby cell, or is it background?",
the answer depends on this rate. The background acts as a constant, location-independent
competitor. A spot only gets assigned to a cell if the cell explains it **better than
pure noise would**.

## Per gene, not one global number

The original pciSeq used a single noise level for every gene. This version learns a
**separate background rate for each gene**, because genes genuinely differ: some are
noisy, some are clean. Letting each gene have its own rate means a noisy gene's spots
are judged more sceptically, while a clean gene's spots are trusted more readily.

The rate is updated every round from the spots that ended up labelled as background:

$$
\text{background rate of gene } g \;\approx\;
\frac{\text{noise spots of gene } g}{\text{area of the tissue}}
$$

This is the first appearance of a shape you will see again and again in pciSeq: an
estimate that is **an observed amount divided by the size of the thing it is spread
over**. Keep an eye out for it, because the next block is built almost entirely out of
ratios like this.

## What feeds in and what comes out

- **Feeds in:** which spots were labelled as background on the previous round.
- **Comes out:** a per-gene background rate, used by
  [block 4](spots-to-cells.md) as the "noise" option each spot is compared against.