---
title: Working with the image
description: Tiling the background image for the viewer, reading it back, and adding the image and the cell boundaries to the SpatialData store.
---

# Working with the image

[`fit`](./reference#fit) never sees the image. It works from the spots and the
segmentation alone. The functions on this page handle the image separately: they tile it
for the viewer, read it back, and add it to the SpatialData store.

All of them take the image as a numpy array: `(H, W)` for a single plane, `(Z, H, W)` for a
stack, or `(Z, H, W, C)` for a stack with channels.

## The `.mbtiles` file

An `.mbtiles` file is a single SQLite database that holds an image as a pyramid of small
tiles at several zoom levels. A viewer asks only for the tiles on screen at the current
zoom, so an image far too large to load at once can be panned and zoomed.
[pciSeq Viewer](../viewer/overview.md) reads it as its background layer. It is written by
[`stage_image`](#stage-image) and read back by [`read_tiles`](#read-tiles).

One file holds one image. To give the viewer more than one, for example DAPI and GCaMP,
call `stage_image` once for each with a different `name`; the viewer offers them as
channels to switch between, see
[Background channels](https://acycliq.github.io/pciSeq_viewer/docs/loading-data.html#channels).

The format is based on [MBTiles](https://github.com/mapbox/mbtiles-spec) by Mapbox, the
format web maps use, and keeps its two tables, `tiles` and `metadata`. It departs from it in
three ways, so it is not a drop-in MBTiles file and a generic MBTiles reader will not
display it correctly:

- **Planes.** `tiles` has one more column, `plane_id`, so one file holds a whole z-stack.
  MBTiles has no such column and holds a single map.
- **Pixels, not geography.** MBTiles tiles are geographic, in the Spherical Mercator
  projection. Here they are in the pixel space of the image and there is no projection, so
  the geographic metadata of MBTiles (`bounds`, `center`) is absent.
- **Rows count from the top.** MBTiles follows the TMS scheme, where `tile_row` 0 is the
  bottom row. Here `tile_row` 0 is the top row, as in the image.

### Tables

`tiles`, one row per tile, with a unique index on the first four columns:

| Column | Type | Meaning |
| --- | --- | --- |
| `plane_id` | integer | The z-plane. 0 for a 2D image. |
| `zoom_level` | integer | 0 is the whole plane in one tile. Each level doubles the resolution. |
| `tile_column` | integer | Counted from the left. |
| `tile_row` | integer | Counted from the top. |
| `tile_data` | blob | The tile, a 256 by 256 pixel JPEG. |

A tile that is blank is not stored, so a plane can hold fewer tiles than its grid.

`metadata`, one row per key, both columns text:

| `name` | `value` |
| --- | --- |
| `format` | `jpg` |
| `minzoom`, `maxzoom` | The lowest and the deepest zoom level. |
| `width`, `height` | Size of the original image in pixels. |
| `plane_count`, `planes` | Number of planes, and their ids as a comma-separated list. |
| `tile_count` | Number of tiles in the file. |
| `created` | When the file was written, UTC. |
| `name`, `description` | As given to `stage_image`. Absent when not given. |
| `tint` | Optional `#RRGGBB` colour the viewer tints the layer with. |

### From a pixel to a tile

Before tiling, each plane is resized so that its longest side is $256 \cdot 2^{\texttt{maxzoom}}$
pixels. At zoom level $z$ a pixel $(x, y)$ of the original image therefore falls in the tile

$$
\texttt{tile\_column} = \Big\lfloor \frac{x\, s_z}{256} \Big\rfloor , \qquad
\texttt{tile\_row} = \Big\lfloor \frac{y\, s_z}{256} \Big\rfloor , \qquad
s_z = \frac{256 \cdot 2^{z}}{\max(\texttt{width}, \texttt{height})} .
$$

Being plain SQLite, the file can be opened directly:

```python
import io, sqlite3
from PIL import Image

con = sqlite3.connect('out/run1/dapi.mbtiles')
meta = dict(con.execute('select name, value from metadata'))

blob, = con.execute('select tile_data from tiles where plane_id = ? and zoom_level = ? '
                    'and tile_column = ? and tile_row = ?', (34, 8, 100, 120)).fetchone()
tile = Image.open(io.BytesIO(blob))        # 256 by 256
```

[`read_tiles`](#read-tiles) does this lookup for a whole region and stitches the tiles.

## `stage_image`

Turns an image or a z-stack into one `.mbtiles` file. Returns the path of the file.

```python
import pciSeq
from skimage import io

dapi = io.imread('dapi.tif')                # (Z, H, W)
path = pciSeq.stage_image(dapi, out_dir='out/run1', name='dapi')
```

`zoom_levels` is the deepest level of the pyramid. Levels 0 to `zoom_levels` are written,
and level $z$ is $256 \cdot 2^z$ pixels on its longer side. By default it is the first
level that holds the image at full resolution: 5 (8192 pixels) for an image 6408 pixels
wide, 6 for one 10000 pixels wide. A deeper level adds no detail, only file size, and
pciSeq Viewer enlarges the deepest level when zooming in past it. See
[`stage_image`](./reference#stage-image).

## `tile_maker`

Makes the same pyramid as plain image files on disk, one folder per plane, without packing
it into an `.mbtiles` file. `stage_image` calls it; it is also useful on its own for any
slippy-map viewer that reads a folder of tiles. The tiling is done by libvips.

```python
info = pciSeq.tile_maker(dapi, out_dir='out/tiles')
info['num_planes'], info['original_dims']
```

The output folder is deleted and recreated if it exists. See
[`tile_maker`](./reference#tile-maker).

## `read_tiles`

Reads a region of a plane back out of an `.mbtiles` file, as an RGB image. `bbox` is in the
pixel coordinates of the original image, the same coordinates as `cellData` and the spots.

```python
im, scale = pciSeq.read_tiles('out/run1/dapi.mbtiles', plane=57,
                              bbox=(5393, 702, 5603, 842), width=1200)
```

`scale` is the zoom factor of the returned image. A point `(x, y)` of the original image is
at `((x - x0) * scale, (y - y0) * scale)` in it, which is what is needed to draw cells or
spots on top. The tiles are JPEG and every level was produced by resizing, so the result is
a close visual copy, not the raw pixels: use it for figures and for checking a segmentation
against the image, not for measurements. `read_tiles` does not need libvips. See
[`read_tiles`](./reference#read-tiles).

## `add_image`

Adds the image to an existing [SpatialData store](./spatialdata-store.md), as one more
element next to the spots, the labels and the tables. Nothing else in the store is touched.

```python
store = 'out/run1/pciSeq/data/spatialdata.zarr'
pciSeq.add_image(store, dapi, name='background')
```

The voxel size is read from the store's own provenance, so the image lines up with the
other elements in microns. See [`add_image`](./reference#add-image).

## `add_boundaries`

Adds the cell outlines to the store as shapes. It needs no data: the outlines are traced
from the segmentation the store already holds.

```python
pciSeq.add_boundaries(store)
```

Shapes in SpatialData are 2D, so a 3D run gets one shapes element per plane,
`cell_boundaries_plane_000` and so on, each indexed by the cell label. See
[`add_boundaries`](./reference#add-boundaries).
