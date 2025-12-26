import shutil
import os
import tempfile
import pyvips
import logging
import numpy as np

from .mbtiles import disk_to_mbtiles

stage_image_logger = logging.getLogger(__name__)

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
        stage_image_logger.info('Moving to the next row: %d/%d '% (j, tiles_down-1) )
        y_top_left = j * tile_size
        for i in range(tiles_across):
            x_top_left = i * tile_size
            tile = image.crop(x_top_left, y_top_left, tile_size, tile_size)
            tile_num = j * tiles_across + i
            fov_id = 'fov_' + str(tile_num)

            out_dir = os.path.join(stage_image_logger.ROOT_DIR, 'fov', fov_id, 'img')
            full_path = os.path.join(out_dir, fov_id +'.tif')
            if not os.path.exists(os.path.dirname(full_path)):
                os.makedirs(os.path.dirname(full_path))
            tile.write_to_file(full_path)
            stage_image_logger.info('tile: %s saved at %s' % (fov_id, full_path) )


def map_image_size(z):
    '''
    returns the image size for each zoom level. Assumes that each map tile is 256x256 pixels
    :param z: 
    :return: 
    '''

    return 256 * 2 ** z


def _process_single_plane(im, zoom_levels, plane_out_dir):
    """Process a single 2D image plane into a tile pyramid.

    Args:
        im: pyvips.Image object
        zoom_levels: number of zoom levels
        plane_out_dir: output directory for this plane's tiles

    Returns:
        pixel_dims: [width, height] after resizing
    """
    # Normalize to 8-bit if not already
    if im.format != 'uchar':
        stage_image_logger.info(f"Converting {im.format} to uchar with normalization")
        mn = im.min()
        mx = im.max()
        
        if mx > mn:
            # Scale to 0-255
            im = (im - mn) * (255.0 / (mx - mn))
        else:
            # Constant image, just offset to 0
            im = im - mn
            
        im = im.cast('uchar')

    dim = map_image_size(zoom_levels)

    # Create output directory
    if not os.path.exists(plane_out_dir):
        os.makedirs(plane_out_dir)

    # The following two lines add an alpha component to rgb which allows for transparency.
    # Is this worth it? It adds quite a bit on the execution time, about x2 increase
    # im = im.colourspace('srgb')
    # im = im.addalpha()

    # Resize to fit the tile pyramid
    factor = dim / max(im.width, im.height)
    im = im.resize(factor)
    stage_image_logger.info('Resized to %d by %d' % (im.width, im.height))
    pixel_dims = [im.width, im.height]

    # Sanity check
    assert max(im.width, im.height) == dim, \
        'Image not scaled properly. Expected %d pixels on longest side' % dim

    # im = im.gravity('south-west', dim, dim) # <---- Uncomment this if the origin is the bottom-left corner

    # Create tile pyramid
    im.dzsave(plane_out_dir, layout='google', suffix='.jpg', background=0)

    return pixel_dims


def tile_maker(img, zoom_levels=8, out_dir=r"./tiles", plane_prefix="plane_"):
    """
    Makes a pyramid of tiles from an image.

    Args:
        img: One of:
            - str: path to a 2D image file (TIFF, PNG, JPEG, etc.)
            - numpy array (H, W): single 2D grayscale image
            - numpy array (Z, H, W): 3D stack of grayscale images
            - numpy array (Z, H, W, C): 3D stack with channels
        zoom_levels: (int) Number of zoom levels to produce. Default is 8.
        out_dir: (str) Output folder for the tile pyramid. Will be deleted and recreated if exists.
        plane_prefix: (str) Prefix for plane subdirectories when processing 3D images.
                      Default is "plane_" resulting in "plane_0", "plane_1", etc.

    Returns:
        dict with keys:
            - 'original_dims': [width, height] of the original input image
            - 'num_planes': number of planes processed
            - 'zoom_levels': number of zoom levels
    """
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    # Determine original dimensions and number of planes
    if isinstance(img, str):
        img = pyvips.Image.new_from_file(img, access='sequential')
        original_dims = [img.width, img.height]
        num_planes = 1
    elif isinstance(img, np.ndarray):
        if img.ndim == 3 and img.shape[-1] <= 4:
            raise ValueError(
                f"Shape {img.shape} looks like (H, W, C). 2D RGB images are not supported. "
                "Convert to grayscale first, or reshape to (C, H, W) to treat channels as planes."
            )
        original_dims = [img.shape[-2], img.shape[-3]] if img.ndim == 4 else [img.shape[-1], img.shape[-2]]
        num_planes = 1 if img.ndim == 2 else img.shape[0]
    else:
        raise TypeError(f"img must be a file path (str) or numpy array, got {type(img)}")

    stage_image_logger.info('Processing %d plane(s), size: %dx%d' % (num_planes, original_dims[0], original_dims[1]))

    # Process each plane
    for z in range(num_planes):
        if num_planes > 1:
            stage_image_logger.info('Plane %d/%d' % (z + 1, num_planes))

        if isinstance(img, pyvips.Image):
            plane = img
        elif img.ndim == 2:
            plane = _numpy_to_vips(img)
        else:
            plane = _numpy_to_vips(img[z])

        _process_single_plane(plane, zoom_levels, os.path.join(out_dir, f"{plane_prefix}{z}"))

    stage_image_logger.info('Done. Pyramid of tiles saved at: %s' % out_dir)

    return {
        'original_dims': original_dims,
        'num_planes': num_planes,
        'zoom_levels': zoom_levels,
    }


def stage_image(img, out_dir=None, zoom_levels=8, name=None, description=None, plane_prefix="plane_", voxel_size=None):
    """
    Process an image into a viewable format (MBTiles).

    This function:
    1. Creates tile pyramids for all planes
    2. Packages tiles into a single MBTiles file
    3. Cleans up temporary tile files

    Args:
        img: One of:
            - numpy array (H, W): single 2D grayscale image
            - numpy array (Z, H, W): 3D stack of grayscale images
            - numpy array (Z, H, W, C): 3D stack with channels
            - str: path to a 2D image file (legacy support)
        out_dir: (str) Output directory for the .mbtiles file. Default: system temp directory.
        zoom_levels: (int) Number of zoom levels to produce. Default is 8.
        name: (str) Short identifier for the dataset. Optional.
                    Example: "WT94_DAPI"
        description: (str) Detailed description of the dataset. Optional.
                    Example: "DAPI background for WT94 mouse cortex, 84 z-planes at 0.9um spacing"
        plane_prefix: (str) Prefix for plane directories. Default is "plane_".
        voxel_size: (list/tuple) Size of a voxel in microns [x, y, z]. Optional.
                    Example: [0.28, 0.28, 0.7] for 0.28 microns in x/y and 0.7 in z.

    Returns:
        None. Check logs for output location.
    """
    # Determine output directory
    if out_dir is None:
        out_dir = tempfile.gettempdir()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Create temporary directory for tiles (in same location as output)
    tiles_dir = os.path.join(out_dir, "_tiles_temp")
    mbtiles_path = os.path.join(out_dir, "output.mbtiles")

    stage_image_logger.info("Starting image processing...")
    stage_image_logger.info("Output directory: %s" % out_dir)

    try:
        # Step 1: Create tile pyramids
        stage_image_logger.info("Step 1/3: Creating tile pyramids...")
        result = tile_maker(img, zoom_levels=zoom_levels, out_dir=tiles_dir, plane_prefix=plane_prefix)

        # Step 2: Package into MBTiles
        stage_image_logger.info("Step 2/3: Packaging tiles into MBTiles...")
        disk_to_mbtiles(
            tiles_dir,
            mbtiles_path,
            format="jpg",
            batch_size=50000,
            width=result['original_dims'][0],
            height=result['original_dims'][1],
            name=name,
            description=description,
            voxel_size=voxel_size,
        )

        # Step 3: Clean up temporary tiles
        stage_image_logger.info("Step 3/3: Cleaning up temporary files...")
        shutil.rmtree(tiles_dir)

    except Exception as e:
        # Clean up on failure
        if os.path.exists(tiles_dir):
            shutil.rmtree(tiles_dir)
        raise

    stage_image_logger.info("Done! MBTiles file created at: %s" % mbtiles_path)

