import shutil
import os
import re
import tempfile
import pyvips
import logging
import numpy as np
from tqdm.auto import tqdm  # picks widget bars in a notebook, text bars in a terminal

from .mbtiles import disk_to_mbtiles, buffer_to_mbtiles

logger = logging.getLogger(__name__)

# Map numpy dtypes to pyvips format strings
DTYPE_TO_VIPS_FORMAT = {
    'uint8': 'uchar',
    'int8': 'char',
    'uint16': 'ushort',
    'int16': 'short',
    'uint32': 'uint',
    'int32': 'int',
    'float32': 'float',
    'float64': 'double',
}


def _numpy_to_vips(arr):
    """Convert a 2D numpy array to a pyvips Image.

    Args:
        arr: 2D numpy array with shape (H, W) or (H, W, C)

    Returns:
        pyvips.Image
    """
    if arr.ndim == 2:
        height, width = arr.shape
        bands = 1
    elif arr.ndim == 3:
        height, width, bands = arr.shape
    else:
        raise ValueError(f"Expected 2D or 3D array, got {arr.ndim}D")

    # Ensure array is contiguous
    arr = np.ascontiguousarray(arr)

    # Get vips format string
    dtype_name = arr.dtype.name
    if dtype_name not in DTYPE_TO_VIPS_FORMAT:
        raise ValueError(f"Unsupported dtype: {dtype_name}. Supported: {list(DTYPE_TO_VIPS_FORMAT.keys())}")

    vips_format = DTYPE_TO_VIPS_FORMAT[dtype_name]

    return pyvips.Image.new_from_memory(arr.tobytes(), width, height, bands, vips_format)

def split_image(im):
    # DEPRECATED to be removed
    '''
    you can just do:
        im.dzsave('./out', suffix='.tif', skip_blanks=-1, background=0, depth='one', overlap=0, tile_size=2000, layout='google')
    to split the image to smaller squares. However you need to write a couple of line to rename and move the file to the correct
    folders
    :param im:
    :return:
    '''
    im = pyvips.Image.new_from_file(im, access='random')
    tile_size = 2000;

    if im.width % tile_size == 0:
        tiles_across = int(im.width / tile_size)
    else:
        tiles_across = im.width // tile_size + 1


    if im.width % tile_size == 0:
        tiles_down = int(im.height/tile_size)
    else:
        tiles_down = im.height // tile_size + 1

    image = im.gravity('north-west', tiles_across * tile_size, tiles_down * tile_size)

    for j in range(tiles_down):
        logger.info('Moving to the next row: %d/%d '% (j, tiles_down-1) )
        y_top_left = j * tile_size
        for i in range(tiles_across):
            x_top_left = i * tile_size
            tile = image.crop(x_top_left, y_top_left, tile_size, tile_size)
            tile_num = j * tiles_across + i
            fov_id = 'fov_' + str(tile_num)

            out_dir = os.path.join(logger.ROOT_DIR, 'fov', fov_id, 'img')
            full_path = os.path.join(out_dir, fov_id +'.tif')
            if not os.path.exists(os.path.dirname(full_path)):
                os.makedirs(os.path.dirname(full_path))
            tile.write_to_file(full_path)
            logger.info('tile: %s saved at %s' % (fov_id, full_path) )


def map_image_size(z):
    '''
    returns the image size for each zoom level. Assumes that each map tile is 256x256 pixels
    :param z: 
    :return: 
    '''

    return 256 * 2 ** z


def _fit_zoom_level(width, height):
    """The first zoom level that is big enough to hold the image.

    Level z is 256 * 2**z pixels on its longer side, so this is the smallest z
    where that covers the longer side of the image. Any level past it would be the
    same picture enlarged, and the viewer does that by itself when zooming in.
    """
    z = 0
    while map_image_size(z) < max(width, height):
        z += 1
    return z


def _get_img_details(img):
    """Determine image dimensions, number of planes, and convert file paths to pyvips.

    Args:
        img: file path (str), or numpy array (2D, 3D, or 4D)

    Returns:
        (img, original_dims, num_planes) where img may have been converted
        from a file path to a pyvips.Image
    """
    if isinstance(img, str):
        img = pyvips.Image.new_from_file(img, access='sequential')
        return img, [img.width, img.height], 1
    elif isinstance(img, np.ndarray):
        if img.ndim == 3 and img.shape[-1] <= 4:
            raise ValueError(
                f"Shape {img.shape} looks like (H, W, C). 2D RGB images are not supported. "
                "Convert to grayscale first, or reshape to (C, H, W) to treat channels as planes."
            )
        original_dims = [img.shape[-2], img.shape[-3]] if img.ndim == 4 else [img.shape[-1], img.shape[-2]]
        num_planes = 1 if img.ndim == 2 else img.shape[0]
        return img, original_dims, num_planes
    else:
        raise TypeError(f"img must be a file path (str) or numpy array, got {type(img)}")


def _prepare_plane(im, zoom_levels):
    """Normalize to 8-bit and resize a pyvips image to fit the tile pyramid.

    Args:
        im: pyvips.Image object
        zoom_levels: the deepest zoom level (levels 0 to zoom_levels get made)

    Returns:
        im: the prepared pyvips.Image (uchar, resized)
    """
    # Normalize to 8-bit if not already
    if im.format != 'uchar':
        logger.debug(f"Converting {im.format} to uchar with normalization")
        mn = im.min()
        mx = im.max()
        if mx > mn:
            im = (im - mn) * (255.0 / (mx - mn))
        else:
            im = im - mn
        im = im.cast('uchar')

    # Resize to fit the tile pyramid
    dim = map_image_size(zoom_levels)
    factor = dim / max(im.width, im.height)
    im = im.resize(factor)
    logger.debug('Resized to %d by %d' % (im.width, im.height))

    assert max(im.width, im.height) == dim, \
        'Image not scaled properly. Expected %d pixels on longest side' % dim

    return im


def _process_single_plane(im, zoom_levels, plane_out_dir):
    """Process a single 2D image plane into a tile pyramid on disk.

    Args:
        im: pyvips.Image object
        zoom_levels: the deepest zoom level (levels 0 to zoom_levels get made)
        plane_out_dir: output directory for this plane's tiles

    Returns:
        pixel_dims: [width, height] after resizing
    """
    im = _prepare_plane(im, zoom_levels)

    if not os.path.exists(plane_out_dir):
        os.makedirs(plane_out_dir)

    im.dzsave(plane_out_dir, layout='google', suffix='.jpg', background=0)

    return [im.width, im.height]


def tile_maker(img, zoom_levels=None, out_dir=r"./tiles", plane_prefix="plane_", progress_bar=None):
    """
    Makes a pyramid of tiles from an image.

    Args:
        img: One of:
            - str: path to a 2D image file (TIFF, PNG, JPEG, etc.)
            - numpy array (H, W): single 2D grayscale image
            - numpy array (Z, H, W): 3D stack of grayscale images
            - numpy array (Z, H, W, C): 3D stack with channels
        zoom_levels: (int) The deepest zoom level. Levels 0 to zoom_levels are written. If not
            given it is the first level big enough to hold the image, see stage_image.
        out_dir: (str) Output folder for the tile pyramid. Will be deleted and recreated if exists.
        plane_prefix: (str) Prefix for plane subdirectories when processing 3D images.
                      Default is "plane_" resulting in "plane_0", "plane_1", etc.
        progress_bar: (tqdm, optional) if given, ticked once per plane instead of
                      logging a per-plane line. Used by stage_image to drive its bars.

    Returns:
        dict with keys:
            - 'original_dims': [width, height] of the original input image
            - 'num_planes': number of planes processed
            - 'zoom_levels': the deepest zoom level that was written
    """
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    img, original_dims, num_planes = _get_img_details(img)
    if zoom_levels is None:
        zoom_levels = _fit_zoom_level(*original_dims)

    logger.info('Processing %d plane(s), size: %dx%d' % (num_planes, original_dims[0], original_dims[1]))

    # Process each plane
    for z in range(num_planes):
        if progress_bar is None and num_planes > 1:
            logger.info('Plane %d/%d' % (z + 1, num_planes))

        if isinstance(img, pyvips.Image):
            plane = img
        elif img.ndim == 2:
            plane = _numpy_to_vips(img)
        else:
            plane = _numpy_to_vips(img[z])

        _process_single_plane(plane, zoom_levels, os.path.join(out_dir, f"{plane_prefix}{z}"))

        if progress_bar is not None:
            progress_bar.update(1)

    logger.info('Done. Pyramid of tiles saved at: %s' % out_dir)

    return {
        'original_dims': original_dims,
        'num_planes': num_planes,
        'zoom_levels': zoom_levels,
    }


def stage_image(img, out_dir=None, zoom_levels=None, name=None, description=None, plane_prefix="plane_",
                use_buffer=True, tint=None, progress=True):
    """
    Turn an image (or z-stack) into an MBTiles file the viewer can read.

    This builds the tiled, multi-resolution background that pciSeq Viewer uses
    as its slippy-map base layer. It takes a single 2D image or a 3D stack and
    writes one `.mbtiles` file.

    Parameters
    ----------
    img : np.ndarray or str
        The image to tile: a 2D array (H, W), a 3D stack (Z, H, W), a 3D stack
        with channels (Z, H, W, C), or a path to a 2D image file (legacy).
    out_dir : str, optional
        Directory for the `.mbtiles` file. Defaults to the system temp directory.
    zoom_levels : int, optional
        The deepest zoom level. Levels 0 to `zoom_levels` are written, and level z
        is 256 * 2**z pixels on its longer side. If not given, it is the first level
        that holds the image at full resolution, for example 5 (8192 pixels) for an
        image 6408 pixels wide. A deeper level adds no detail, only file size; the
        viewer enlarges the deepest level when zooming in past it.
    name : str, optional
        Short identifier for the dataset. Also used as the output filename, e.g.
        `name="S10_gcamp_10"` writes `S10_gcamp_10.mbtiles`. If empty, the file
        is named `output.mbtiles`.
    description : str, optional
        Longer description of the dataset.
    plane_prefix : str, optional
        Prefix for the per-plane names. Default is "plane_".
    use_buffer : bool, optional
        If True (the default) the tiles are built in memory and inserted straight
        into the MBTiles database. If False they are written to disk first, which
        uses less memory but more disk I/O.
    tint : str, optional
        Hex colour like "#00FF00" the viewer uses to tint this grayscale layer.
        If omitted, the layer is shown in plain grayscale.
    progress : bool, optional
        If True (the default) show two per-plane tqdm bars, one for tiling and
        one for the db writing. Set to False for headless/quiet runs.

    Returns
    -------
    str
        Path to the created `.mbtiles` file.

    Notes
    -----
    Requires libvips. If it is not installed, `pciSeq.stage_image()` raises an
    ImportError that explains how to install it.
    """
    # Determine output directory
    if out_dir is None:
        out_dir = tempfile.gettempdir()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # name the output file after `name` so it is self-identifying and the viewer
    # can use the filename stem as the channel id. Non filename-safe characters
    # are replaced with '_'. If `name` is empty, fall back to output.mbtiles.
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()) if name and name.strip() else ""
    mbtiles_path = os.path.join(out_dir, f"{safe_name}.mbtiles" if safe_name else "output.mbtiles")

    logger.info("Starting image processing...")
    logger.info("Output directory: %s" % out_dir)

    if zoom_levels is None:
        # stop at the level that holds the image, deeper ones are only enlargements
        zoom_levels = _fit_zoom_level(*_get_img_details(img)[1])
        logger.info("Deepest zoom level: %d" % zoom_levels)

    if use_buffer:
        _stage_image_buffer(img, mbtiles_path, zoom_levels, name, description, plane_prefix, tint, progress)
    else:
        _stage_image_disk(img, mbtiles_path, zoom_levels, name, description, plane_prefix, tint, progress)

    logger.info("Done! MBTiles file created at: %s" % mbtiles_path)

    return mbtiles_path


def _plane_buffer_generator(img, num_planes, zoom_levels, plane_prefix, tile_bar=None):
    """Yield one dzsave_buffer per plane. O(1) memory, only one plane's tiles in memory at a time.

    If tile_bar is given it is ticked once per plane, right after that plane's
    tiles are built, so the caller's "making tiles" bar stays in step with the work.
    """
    for z in range(num_planes):
        if tile_bar is None and num_planes > 1:
            logger.info('Plane %d/%d' % (z + 1, num_planes))

        if isinstance(img, pyvips.Image):
            plane = img
        elif img.ndim == 2:
            plane = _numpy_to_vips(img)
        else:
            plane = _numpy_to_vips(img[z])

        plane = _prepare_plane(plane, zoom_levels)
        buf = plane.dzsave_buffer(basename=f'{plane_prefix}{z}', layout='google', suffix='.jpg', background=0)
        if tile_bar is not None:
            tile_bar.update(1)
        yield buf


def _stage_image_buffer(img, mbtiles_path, zoom_levels, name, description, plane_prefix, tint=None, progress=True):
    """In-memory path: tiles never touch disk. The SQLite db is built in a
    local temp directory and copied to mbtiles_path at the end, so that
    SQLite never opens the db on a network filesystem (NFS/SMB locking is
    unreliable and can corrupt the file)."""
    img, original_dims, num_planes = _get_img_details(img)

    logger.info('Processing %d plane(s), size: %dx%d' % (num_planes, original_dims[0], original_dims[1]))
    logger.info("Creating tile pyramids and packaging into MBTiles...")

    # Two per-plane bars. Because bufs is a lazy generator, tiling and db writing
    # happen one plane at a time, so the two bars advance almost together (the db
    # bar trails the tile bar by ~1 plane).
    with tempfile.TemporaryDirectory() as tmpdir, \
            tqdm(total=num_planes, desc="Making tiles", unit="plane", position=0, disable=not progress) as tile_bar, \
            tqdm(total=num_planes, desc="Writing DB  ", unit="plane", position=1, disable=not progress) as db_bar:
        local_mbtiles = os.path.join(tmpdir, "output.mbtiles")
        bufs = _plane_buffer_generator(img, num_planes, zoom_levels, plane_prefix, tile_bar=tile_bar)
        buffer_to_mbtiles(
            bufs,
            local_mbtiles,
            format="jpg",
            batch_size=50000,
            width=original_dims[0],
            height=original_dims[1],
            name=name,
            description=description,
            tint=tint,
            plane_bar=db_bar,
        )
        shutil.copy2(local_mbtiles, mbtiles_path)


def _stage_image_disk(img, mbtiles_path, zoom_levels, name, description, plane_prefix, tint=None, progress=True):
    """Disk-based path: tiles AND db are built in a local temp directory,
    then the finished db is copied to mbtiles_path. Local-first avoids
    SQLite file-locking issues on network filesystems (NFS/SMB), and the
    temp directory is auto-cleaned on exit (including on exception)."""
    # plane count up front so we can size both bars. On this path the tile bar
    # fills completely first, then the db bar, since all tiles are written to
    # disk before any get read back into the db.
    num_planes = _get_img_details(img)[2]

    with tempfile.TemporaryDirectory() as tmpdir, \
            tqdm(total=num_planes, desc="Making tiles", unit="plane", position=0, disable=not progress) as tile_bar, \
            tqdm(total=num_planes, desc="Writing DB  ", unit="plane", position=1, disable=not progress) as db_bar:
        tiles_dir = os.path.join(tmpdir, "tiles")
        local_mbtiles = os.path.join(tmpdir, "output.mbtiles")

        result = tile_maker(img, zoom_levels=zoom_levels, out_dir=tiles_dir,
                            plane_prefix=plane_prefix, progress_bar=tile_bar)

        disk_to_mbtiles(
            tiles_dir,
            local_mbtiles,
            format="jpg",
            batch_size=50000,
            width=result['original_dims'][0],
            height=result['original_dims'][1],
            name=name,
            description=description,
            tint=tint,
            plane_bar=db_bar,
        )

        shutil.copy2(local_mbtiles, mbtiles_path)

