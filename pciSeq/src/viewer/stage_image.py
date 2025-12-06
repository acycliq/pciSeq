import shutil
import os
import pyvips
import logging
import numpy as np

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
            - numpy array (H, W) or (H, W, C): single 2D image
            - numpy array (Z, H, W) or (Z, H, W, C): 3D stack of images (multiple planes)
        zoom_levels: (int) Number of zoom levels to produce. Default is 8.
        out_dir: (str) Output folder for the tile pyramid. Will be deleted and recreated if exists.
        plane_prefix: (str) Prefix for plane subdirectories when processing 3D images.
                      Default is "plane_" resulting in "plane_0", "plane_1", etc.

    Returns:
        dict with keys:
            - 'pixel_dims': [width, height] of the resized image
            - 'num_planes': number of planes processed
            - 'zoom_levels': number of zoom levels
    """
    # Remove output dir if exists, then create fresh
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    # Determine input type and process accordingly
    if isinstance(img, str):
        # File path - load as single 2D image
        stage_image_logger.info('Loading image from file: %s' % img)
        im = pyvips.Image.new_from_file(img, access='sequential')
        plane_out_dir = os.path.join(out_dir, f"{plane_prefix}0")
        stage_image_logger.info('Processing single image...')
        pixel_dims = _process_single_plane(im, zoom_levels, plane_out_dir)
        num_planes = 1

    elif isinstance(img, np.ndarray):
        # Numpy array - check dimensionality
        if img.ndim == 2 or (img.ndim == 3 and img.shape[2] <= 4):
            # 2D image: (H, W) or (H, W, C) where C is channels (1-4)
            stage_image_logger.info('Processing 2D numpy array with shape %s' % (img.shape,))
            im = _numpy_to_vips(img)
            plane_out_dir = os.path.join(out_dir, f"{plane_prefix}0")
            pixel_dims = _process_single_plane(im, zoom_levels, plane_out_dir)
            num_planes = 1

        elif img.ndim == 3 or img.ndim == 4:
            # 3D stack: (Z, H, W) or (Z, H, W, C)
            num_planes = img.shape[0]
            stage_image_logger.info('Processing 3D numpy array with %d planes, shape %s' % (num_planes, img.shape))

            pixel_dims = None
            for z in range(num_planes):
                stage_image_logger.info('Processing plane %d/%d' % (z + 1, num_planes))
                plane_data = img[z]
                im = _numpy_to_vips(plane_data)
                plane_out_dir = os.path.join(out_dir, f"{plane_prefix}{z}")
                pixel_dims = _process_single_plane(im, zoom_levels, plane_out_dir)

        else:
            raise ValueError(f"Unsupported numpy array dimensions: {img.ndim}")

    else:
        raise TypeError(f"img must be a file path (str) or numpy array, got {type(img)}")

    stage_image_logger.info('Done. Pyramid of tiles saved at: %s' % out_dir)

    return {
        'pixel_dims': pixel_dims,
        'num_planes': num_planes,
        'zoom_levels': zoom_levels,
    }



