# tiling

Turns a background image into map tiles for the viewer. This is image prep
only, it is not part of the cell-typing algorithm and nothing in the core
pipeline imports it. It builds the tiles a viewer reads.

## Files

- `stage_image.py` - the public entry point. Two functions:
  - `tile_maker(img, ...)` builds a zoom-level tile pyramid from an image and
    writes the tiles to a folder on disk.
  - `stage_image(img, ...)` does the same but packs the result into a single
    `.mbtiles` file (an SQLite database of tiles) and returns its path. The file
    is named after `name` (e.g. `name="S10_gcamp_10"` -> `S10_gcamp_10.mbtiles`),
    or `output.mbtiles` when `name` is empty.
- `mbtiles.py` - the helper that `stage_image` uses to write tiles into the
  `.mbtiles` SQLite format (`disk_to_mbtiles`, `buffer_to_mbtiles`). It is only
  imported by `stage_image.py`.

## How it is used

`tile_maker` and `stage_image` are re-exported at the top level, so callers use
`pciSeq.tile_maker(...)` / `pciSeq.stage_image(...)`. Both need libvips
(via `pyvips`). If libvips is missing, these become stubs that just warn, so
importing pciSeq still works.

## Background tint

`stage_image(..., tint="#RRGGBB")` writes a `tint` key into the `.mbtiles`
`metadata` table. The tiles themselves stay grayscale; the hex is just a hint the
viewer reads back to colour that background channel (e.g. green for a GCaMP
channel). If `tint` is omitted or not a valid `#RRGGBB` string, no key is written
and the viewer renders the layer grayscale. The tint is a property of the
channel, not the pixels, so you pass it in rather than deriving it from the
(single-channel, grayscale) image.

