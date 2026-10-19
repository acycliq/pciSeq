"""Read a region back out of an .mbtiles pyramid.

stage_image goes one way, image to tiles. This goes the other way: give it a box in
image pixels and it hands back that piece of the picture. It only reads the tiles that
touch the box, so cropping a cell out of a 10GB file is a couple of tile reads, not a
full stitch.

Needs nothing but sqlite3 and pillow, so it works even where libvips is missing.
"""

import io
import math
import sqlite3
from typing import Optional, Tuple

from PIL import Image

TILE = 256


def _metadata(con) -> dict:
    return {k: v for k, v in con.execute("select name, value from metadata")}


def _level_scale(zoom: int, width: int, height: int) -> float:
    """Tile pixels per image pixel at this zoom level.

    stage_image resizes every plane so its longest side is TILE * 2 ** zoom_levels,
    which is where this comes from.
    """
    return TILE * 2 ** zoom / max(width, height)


def read_tiles(mbtiles: str,
               plane: int = 0,
               bbox: Optional[Tuple[float, float, float, float]] = None,
               width: Optional[int] = None,
               zoom: Optional[int] = None) -> Tuple[Image.Image, float]:
    """Read part of a plane out of an mbtiles pyramid.

    Parameters
    ----------
    mbtiles : str
        Path to the .mbtiles file, as written by :func:`pciSeq.stage_image`.
    plane : int, default 0
        The z-plane to read. For a 2D image there is only plane 0.
    bbox : tuple of float, optional
        The region to read, ``(x0, y0, x1, y1)`` in the pixel coordinates of the
        original image, the same coordinates as `cellData` and the spots. The box is
        clamped to the image. The default is the whole plane.
    width : int, optional
        Width of the returned image in pixels. The pyramid level is picked to cover it
        and the result is resampled to exactly this width. The default returns the
        region at its original scale, one pixel per image pixel.
    zoom : int, optional
        Read this pyramid level instead of choosing one. Mostly useful for inspecting
        the file itself.

    Returns
    -------
    image : PIL.Image.Image
        The region, in RGB.
    scale : float
        Zoom factor, returned width over box width. A point ``(x, y)`` of the original
        image is at ``((x - x0) * scale, (y - y0) * scale)`` in the returned one.

    Notes
    -----
    This is not a lossless recovery of the original data. The tiles are JPEG and every
    level was produced by resizing, so the crop is a close visual copy, not the raw
    pixels of the image the pyramid was built from. Use it for figures and for checking
    a segmentation against the image, not for measurements.

    `width` sets the size of the result and therefore how much detail is read. The
    pyramid level is chosen as the cheapest one that covers the request, and the result
    is resampled to exactly `width`. Levels only go as high as the pyramid does, so a
    `width` beyond the top level is allowed and enlarges the result, which is
    magnification, not detail.

    Examples
    --------
    A cell and its surroundings, one pixel per image pixel:

    >>> im, scale = read_tiles('dapi.mbtiles', plane=57, bbox=(5393, 702, 5603, 842))
    >>> im.size, scale
    ((210, 140), 1.0)

    The same region as a 1200 pixel wide panel. The box is 210 across, so `scale` is
    5.71 and the panel is magnified: the level read is the top of the pyramid, which
    renders that box at 10.2 times magnification, and it is resampled down to 5.71.

    >>> im, scale = read_tiles('dapi.mbtiles', plane=57, bbox=(5393, 702, 5603, 842),
    ...                        width=1200)

    A whole 6408 by 4382 plane, trimmed to 3:2 and shown 1200 wide. Here `scale` is
    1200 / 6408 = 0.187, and a cell at (5498.2, 772.5) lands at
    ``((5498.2 - 0) * 0.187, (772.5 - 55) * 0.187)``, so about (1029, 134) in the
    returned image. `bbox=None` would give the untrimmed plane.

    >>> im, scale = read_tiles('dapi.mbtiles', plane=57, bbox=(0, 55, 6408, 4327),
    ...                        width=1200)
    """
    con = sqlite3.connect(f"file:{mbtiles}?mode=ro", uri=True)
    try:
        meta = _metadata(con)
        img_w, img_h = int(meta["width"]), int(meta["height"])
        maxzoom = int(meta["maxzoom"])

        x0, y0, x1, y1 = (0.0, 0.0, float(img_w), float(img_h)) if bbox is None else map(float, bbox)
        if x1 <= x0 or y1 <= y0:
            raise ValueError(f"bbox must be (x0, y0, x1, y1) with x1 > x0 and y1 > y0, got {bbox}")
        x0, x1 = max(0.0, x0), min(float(img_w), x1)
        y0, y1 = max(0.0, y0), min(float(img_h), y1)
        box_w, box_h = x1 - x0, y1 - y0

        # how many returned pixels per image pixel the caller asked for
        want = width / box_w if width else 1.0

        if zoom is None:
            # the cheapest level that still has the detail asked for
            zoom = next((z for z in range(maxzoom + 1) if _level_scale(z, img_w, img_h) >= want),
                        maxzoom)
        s = _level_scale(zoom, img_w, img_h)

        # the box in the tile pixels of this level, and the tiles it touches
        tx0, ty0 = math.floor(x0 * s / TILE), math.floor(y0 * s / TILE)
        tx1, ty1 = math.floor((x1 * s - 1e-6) / TILE), math.floor((y1 * s - 1e-6) / TILE)
        canvas = Image.new("RGB", ((tx1 - tx0 + 1) * TILE, (ty1 - ty0 + 1) * TILE))
        rows = con.execute(
            "select tile_column, tile_row, tile_data from tiles where plane_id = ? "
            "and zoom_level = ? and tile_column between ? and ? and tile_row between ? and ?",
            (plane, zoom, tx0, tx1, ty0, ty1))
        found = 0
        for tx, ty, blob in rows:
            canvas.paste(Image.open(io.BytesIO(blob)).convert("RGB"),
                         ((tx - tx0) * TILE, (ty - ty0) * TILE))
            found += 1
        if found == 0:
            raise ValueError(f"no tiles for plane {plane} at zoom {zoom} in {mbtiles}")

        im = canvas.crop((round(x0 * s - tx0 * TILE), round(y0 * s - ty0 * TILE),
                          round(x1 * s - tx0 * TILE), round(y1 * s - ty0 * TILE)))

        out_w, out_h = max(1, round(box_w * want)), max(1, round(box_h * want))
        if (im.width, im.height) != (out_w, out_h):
            im = im.resize((out_w, out_h), Image.LANCZOS)
        return im, im.width / box_w
    finally:
        con.close()
