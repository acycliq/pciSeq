---
title: pciSeq Viewer
description: "The desktop viewer for a finished run: the 2D map, the 3D voxel view, and other features."
---

# pciSeq Viewer

pciSeq Viewer is a desktop application that reads the output of a finished run. It draws
the spots and the cells over the background image, one plane at a time, and reads the same
diagnostics as [`check_cell`](../api/reference.md#check-cell) and
[`check_spot`](../api/reference.md#check-spot), so a call can be inspected by Ctrl+Click on it.

It is a separate project, with its own
[documentation](https://acycliq.github.io/pciSeq_viewer/) and installers for Windows,
macOS and Linux on its
[releases page](https://github.com/acycliq/pciSeq_viewer/releases/latest). A web version
that requires no installation is at
[web viewer](https://acycliq.github.io/pciSeq_viewer/web_viewer/).

It is not the [live viewer](../api/live-viewer.md). The live viewer runs in a browser while
`fit` is running and shows the class assignments as they are updated. pciSeq Viewer is
opened after the run, on the saved output.

## Input files

It reads two inputs:

- **`viewer_data/`**, written by a run with `save_data` on, under `<output_path>/pciSeq/data`, see
  [working with results](../api/working-with-results.md#saved-files). It holds the spots and
  the cell boundaries as feather files, one per plane, the cells as feather files of up to
  100,000 rows each, and `diagnostics/diagnostics.db`.
- **An `.mbtiles` file**, the tiled background image, written by
  [`stage_image`](../api/working-with-the-image.md#stage-image). Without it the viewer asks for
  the image dimensions and draws everything on a blank background.

The viewer is opened on the `viewer_data` folder and locates the other files from it.

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

The main window also provides controls to show or hide genes and classes, filter the spots by score or intensity, change
the opacity of each layer, switch background channel when the image has more than one,
draw a region and export it, and open the class and misread charts.

## The 3D voxel view

Selecting a rectangle on the map opens a second window that renders that region of tissue in
3D: the cells, their boundaries and the spots become voxels in a scene that can be
orbited, rather than a stack of planes. It is restricted to the selected region, which
limits the rendering cost on large datasets.

## The diagnostics

With `diagnostics/diagnostics.db` connected, Ctrl+Click a cell to open the cell diagnostics
and Ctrl+Click a spot to open the spot diagnostics. They show the same terms as
[`check_cell`](../api/reference.md#check-cell) and [`check_spot`](../api/reference.md#check-spot), the
first comparing the assigned class against another class selected by the user, the second the
score of every candidate cell of the spot. The two pages under
[explaining the calls](../explaining-the-calls/overview.md) describe those terms in
detail.
