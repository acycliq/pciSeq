"""
Functions to rotate an image and then make a pyramid of tiles
"""

import numpy as np
import pandas as pd
import os
import pyvips
import pciSeq.src.preprocess.vizgen.config as config
import shutil
import logging
from pathlib import Path

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def rotate_image(merfish_id, slice_id, region_id):
    """
    rotates an image.
    img_in: path to the image to be rotated
    img_out: path to save the rotared image to
    deg: degrees to rotate the image by (clockwise)
    """
    cfg = config.get_config(merfish_id=merfish_id, slice_id=slice_id, region_id=region_id)
    deg = cfg['rotation'][0]
    logger.info('Rotating %s\%s by %d degrees' % (slice_id, region_id, deg))

    img_in = os.path.join('Z:\\', merfish_id, slice_id, region_id, 'images', 'mosaic_DAPI_z3.tif')
    img_out = os.path.join('Z:\\', 'Dimitris folder', 'rotated_dapi', slice_id, region_id, 'images', 'mosaic_DAPI_z3.tif')

    if not os.path.exists(img_out):
        _dir = Path(img_out).parent.absolute()
        os.makedirs(_dir)

    x = pyvips.Image.new_from_file(img_in)
    x = x.rotate(deg, interpolate=pyvips.Interpolate.new("nearest"))
    x.write_to_file(img_out, compression="jpeg", tile=True)

    return img_in, img_out


def map_image_size(z):
    '''
    return the image size for each zoom level. Assumes that each map tile is 256x255
    :param z:
    :return:
    '''

    return 256 * 2 ** z


def tile_maker(z_depth, out_dir, img_path):
    """
    makes a pyramid of map tiles.
    :param z_depth: the maximum zoom level
    :param out_dir: the path to the folder where the map tiles will be saved
    :param img_path: the full path of the image to be tiled
    :return:
    """
    logger.info('Making %d-deep pyramid of map tiles at %s' % (z_depth, out_dir))

    dim = map_image_size(z_depth)
    # remove the dir if it exists
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    # now make a fresh one
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    im = pyvips.Image.new_from_file(img_path, access='sequential')

    # The following two lines add an alpha component to rgb which allows for transparency.
    # Is this worth it? It adds quite a bit on the execution time, about x2 increase
    # im = im.colourspace('srgb')
    # im = im.addalpha()

    logger.info('Resizing image: %s' % img_path)
    factor = dim / max(im.width, im.height)
    im = im.resize(factor)
    logger.info('Done! Image is now %d by %d' % (im.width, im.height))
    pixel_dims = [im.width, im.height]

    # sanity check
    assert max(im.width, im.height) == dim, 'Something went wrong. Image isnt scaled up properly. ' \
                                            'It should be %d pixels in its longest side' % dim

    # im = im.gravity('south-west', dim, dim) # <---- Uncomment this if the origin is the bottomleft corner

    # now you can create a fresh one and populate it with tiles
    logger.info('Started doing the image tiles ')
    im.dzsave(out_dir, layout='google', suffix='.jpg', background=0, skip_blanks=0)
    logger.info('Done. Pyramid of tiles saved at: %s' % out_dir)

    return pixel_dims


def run(merfish_id, slice_id, region_id):
    img_in, img_out = rotate_image(merfish_id, slice_id, region_id)
    _dir = os.path.join(os.path.dirname(img_out), '262144px')
    # tile_maker(10, _dir, img_out)



if __name__ == "__main__":
    merfish_id = 'MERFISH_F_E'
    slice_ids = [
        "VS6_MsBrain_B3_VS6library_V10_LH_02-07-21",
        # "MsBrain_ZM1_VS6_JH_V11_05-16-2021",
        # "MsBrain_ZM2_VS6_JH_V11_05-15-2021",
        # "MsBrain_ZM3_VS6_JH_V11_05-17-2021",
        # "MsBrain_ZM4_VS6_JH_V11_05-11-2021",
        # "MsBrain_ZM5.1_VS6_JH_V11_05-12-2021",
        # "MsBrain_ZM5.2_VS6_JH_V6_05-13-2021",
        # "MsBrain_ZM6.1_VS6_V6_JH_05-11-2021",
        # "MsBrain_ZM7.1_VS6_V6_JH_05-12-2021",
        # "MsBrain_ZM7.2_VS6_JH_V11_05-13-2021"
        ]
    region_ids = ['region_0', 'region_1']

    for slice_id in slice_ids:
        for region_id in region_ids:
            run(merfish_id, slice_id, region_id)

    logger.info('Done!')