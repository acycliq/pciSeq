---
title: pciSeq Viewer
description: "The desktop viewer for a finished run: the 2D map, the 3D voxel view, and other features."
---

# pciSeq Viewer

pciSeq Viewer is a desktop application that reads the output of a finished run. It draws
the spots and the cells over the background image, one plane at a time, and reads the same
diagnostics as [`check_cell`](../api/reference.md#check-cell) and
[`check_spot`](../api/reference.md#check-spot), so a call can be inspected by clicking on it.

It is a separate project, with its own
[documentation](https://acycliq.github.io/pciSeq_viewer/) and installers for Windows,
macOS and Linux on its
[releases page](https://github.com/acycliq/pciSeq_viewer/releases/latest). A web version
that needs no install is at
[web viewer](https://acycliq.github.io/pciSeq_viewer/web_viewer/).

It is not the [live viewer](../api/live-viewer.md). The live viewer runs in a browser while
`fit` is still going and shows the classes as they settle. This one is opened afterwards,
on the saved output.

## What it reads

Two things, both written by a run with `save_data` on:

- **`viewer_data/`**, under `output_path`, see
  [working with results](../api/working-with-results.md#saved-files). It holds the spots, the
  cells and the cell boundaries as feather files, one per plane, and
  `diagnostics/diagnostics.db`.
- **An `.mbtiles` file**, the tiled background image, written by
  [`stage_image`](../api/working-with-the-image.md#stage-image). Without it the viewer asks for
  the image dimensions and draws everything on a blank background.

Point the viewer at the `viewer_data` folder and it picks up the rest.

The folder can also be given on the command line with `--data`. The viewer then opens the
run directly, without the welcome screen:

```bash
pciSeq_viewer --data=/path/to/viewer_data
```

`pciSeq_viewer` stands for the installed executable, whose location depends on the
platform. From a source checkout the command is `npm start -- --data=/path/to/viewer_data`.

## The 2D map

The main window is a map of one plane: the background image underneath, the cell outlines
filled with the colour of the class each cell was assigned to, and the spots as one glyph
per gene. The plane slider moves through the stack.

From there: show or hide genes and classes, filter the spots by score or intensity, change
the opacity of each layer, switch background channel when the image has more than one,
draw a region and export it, and open the class and misread charts.

## The 3D voxel view

Selecting a rectangle on the map opens a second window that builds that piece of tissue in
3D: the cells, their boundaries and the spots become voxels in a scene that can be
orbited, rather than a stack of planes. It works on the selected region only, which is
what keeps it fast on a large dataset.

## The inspectors

With `diagnostics/diagnostics.db` connected, Ctrl-click a cell to open the cell inspector
and Ctrl-click a spot to open the spot inspector. They show the same terms as
[`check_cell`](../api/reference.md#check-cell) and [`check_spot`](../api/reference.md#check-spot), the
first comparing the assigned class against another class of your choice, the second the
score of every candidate cell of the spot. The two pages under
[explaining the calls](../explaining-the-calls/overview.md) walk through those terms in
detail.
