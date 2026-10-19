---
title: Live viewer
description: Watch the cell typing in a browser while pciSeq.fit() runs, and the controls of the page.
---

# Live viewer

The live viewer shows the cells in a browser while [`fit`](./reference#fit) is still
running. Every cell is drawn at its position and coloured by the class it is currently
assigned to, and the picture updates after each iteration. It shows which cells keep
changing class and how quickly a run converges.

The viewer reads the class probabilities after each iteration and does not change them.
A run gives the same result with the viewer on or off.

It is not [pciSeq Viewer](../viewer/overview.md), the desktop application. That one is
opened after a run, on the saved output.

<figure class="diagram">
<!-- the mp4 is hosted on github (issue #1, "docs assets") so it stays out of the repo -->
<video src="https://github.com/user-attachments/assets/ab9867e4-866e-4f81-88e5-359bdf7d7a64"
       poster="/live-viewer/live-viewer-poster.jpg"
       autoplay loop muted playsinline controls width="1200"></video>
<figcaption>The first 14 iterations of a 3D dataset of 25,254 cells, recorded during the run,
with the intervals between iterations removed. Cells are coloured by their current class, the
counts on the right follow the class sizes, and the chart at the bottom right is the
convergence value against its tolerance.</figcaption>
</figure>

## Turning it on

```python
opts = {
    "realtime_viewer": True,
    "realtime_viewer_port": 5001,   # the default
}
cellData, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=scRNAseq, opts=opts)
```

A browser tab opens at `http://127.0.0.1:5001`. If it does not, open that address.
The page says *Waiting for algorithm to start...* until the first iteration is done;
preprocessing runs first, and on a large dataset that takes minutes.

`realtime_viewer_fixed_radius` draws every cell at the given radius. `None`, the default,
draws each cell at the radius derived from its area.

When `fit` returns the viewer shuts down and the page shows it is no longer connected. The
last picture stays, hover and class hiding still work, and anything that needs the
running model, such as the cell diagnostics below, stops.

::: tip Running on a remote machine
The viewer listens on `127.0.0.1` only. To reach it from another machine, forward the
port over ssh and open the address locally:

```bash
ssh -L 5001:127.0.0.1:5001 user@server
```
:::

## The page

**Header.** Whether the page is connected to the run, the current iteration, the
convergence value (see below) and the number of cells. An overlay on the plot also shows
the connection status and the cell count.

**Hover** over a cell for its label, class, position and the probability of that class.

**Pan and zoom** with the mouse.

**Cells that changed class** in the last iteration are enlarged for two seconds.

## Cell classes

The legend lists every class with its colour and how many cells are assigned to it.

- Click a class to hide it, click again to show it. **Show All** and **Hide All** act on
  every class.
- `/` moves the focus to the filter box; typing part of a class name shortens the list.
  `Esc` clears it.
- The handle under the list resizes it.

## Updates

**All Cells** draws every cell. **Updated Only** draws only the cells whose highest class
probability changed by at least the threshold since the previous iteration. The slider
sets the threshold, 1% by default. Two counts show how many cells were updated and how many
changed class. On the first iteration every cell counts as updated.

**Plane** is for 3D data. It adds a slider at the bottom of the picture and draws only the
cells whose centroid falls in the chosen z plane. It starts on the middle plane. The button
is greyed out for 2D data.

## Convergence

The chart plots the convergence value for every iteration: the largest change in any
spot's assignment probability, to a cell or to the background, since the previous iteration. The dashed line is
`CellCallTolerance`, and the run stops once the value drops below it (or at `max_iter`).
Hover over a point to read its value.

## Cell diagnostics

`Ctrl`+click a cell (`Cmd`+click on a Mac) opens the diagnostics drawer. It lists, gene by
gene, the evidence for the cell's class against another class, chosen in **Compare
against**. The default is Zero.

- The two charts show the 10 genes that favour the assigned class the most and the 10
  that favour the compared class the most, measured as the difference in log-likelihood.
- The table has the counts of those genes in this cell next to the mean counts in the
  cells typed as each of the two classes.

Selecting the assigned class in **Compare against** gives an error in the drawer. The top edge of the
drawer resizes it.

## Custom colours

By default every class gets its own colour, spread evenly over the hue wheel, however many
classes there are. `Zero` is always black. **Import Color Scheme** loads a JSON file mapping
class names to colours:

```json
{
  "037 DG Glut": "#5C79CC",
  "038 DG-PIR Ex IMN": "#86FF4D",
  "Zero": "#000000"
}
```

- Names must match the columns of the `scRNAseq` DataFrame exactly.
- A colour is a hex code or a CSS colour name.
- Once a file is loaded, a class it does not name is drawn grey. Those classes are listed
  in the browser console, as are names in the file that are not in the data.
- `Zero` stays black whatever the file says.
- A file loaded before the first iteration is applied once the class names arrive.

::: tip Chrome and Opera on Linux
In the file dialog, select the file and press **Open**. Double-clicking the file does not
always pass it to the page in Chromium based browsers.
:::

## Opening the page during a run

The server keeps the cell positions and the latest iteration, so a tab opened or reloaded
during a run receives the current state immediately. Closing the tab does not affect the run.
