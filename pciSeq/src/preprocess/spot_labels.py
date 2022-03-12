"""
Functions to prepare the data for pciSeq. The label image and spots are parsed and if a spot
lies within the cell boundaries then the corresponding cell id is recorded.
Cell centroids and cell areas are also calculated.
"""

import os
import shutil
import numpy as np
import pandas as pd
import skimage.measure as skmeas
from typing import Tuple
from scipy.sparse import coo_matrix, csr_matrix, save_npz, load_npz
from pciSeq.src.preprocess.cell_borders import extract_borders_par, extract_borders_dip
from PIL import Image, ImageOps, ImageDraw
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger(__name__)


def inside_cell(label_image, spots) -> np.array:
    if isinstance(label_image, coo_matrix):
        label_image = label_image.tocsr()
    elif isinstance(label_image, np.ndarray):
        label_image = csr_matrix(label_image)
    elif isinstance(label_image, csr_matrix):
        pass
    else:
        raise Exception('label_image should be of type "csr_matrix" ')
    m = label_image[spots.y, spots.x]
    out = np.asarray(m)
    return out[0]


def remap_labels(coo):
    """
    Used for debugging/sanity checking only. It resuffles the label_image
    """
    coo_max = coo.data.max()
    _keys = 1 + np.arange(coo_max)
    _vals = _keys.copy()
    np.random.shuffle(_vals)
    d = dict(zip(_keys, _vals))
    new_data = np.array([d[x] for x in coo.data]).astype(np.uint64)
    out = coo_matrix((new_data, (coo.row, coo.col)), shape=coo.shape)
    return out


def filter_data(spots, coo):
    assert np.all([coo[i].shape == coo[i+1].shape for i in range(len(coo)-1)])
    Sx = coo[0].shape[1]
    Sy = coo[0].shape[0]
    Sz = len(coo)
    mask_x = (spots.x >= 0) & (spots.x <= Sx)
    mask_y = (spots.y >= 0) & (spots.y <= Sy)
    mask_z = (spots.z >= 0) & (spots.z <= Sz)
    return spots[mask_x & mask_y & mask_z]


def resize_coo(coo, width, height):
    out = []
    for c in coo:
        _img = np.array(Image.fromarray(c.toarray()).resize((width, height), Image.NEAREST), dtype=np.uint32)
        out.append(coo_matrix(_img))
    return out


def clip_z(coo, spots, z_start, z_end):
    coo_out = coo[z_start: z_end]
    spots_out = spots[(spots.z >= z_start) & (spots.z < z_end)]
    spots_out.z = spots_out.z-z_start
    return coo_out, spots_out


def label_spots(spots, coo):
    """
    assigns a label to the spots. The label is the cell_id within which the spots lies
    If a spot is on the background the label of the spot is zero
    """
    spots_list = []
    for i, d in enumerate(np.unique(spots.z)):
        # z_plane = spots.z == d
        z_plane = int(np.floor(d))  # get the closest z-plane
        spots_z = spots[spots.z == d]  # get all the spots with the same z-coord as d
        my_coo = coo[z_plane]
        inc = inside_cell(my_coo.tocsr(), spots_z)
        spots_list.append(spots_z.assign(label=inc))
    return pd.concat(spots_list)


def stage_data(spots: pd.DataFrame, coo: coo_matrix, cfg) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Reads the spots and the label image that are passed in and calculates which cell (if any) encircles any
    given spot within its boundaries. It also retrieves the coordinates of the cell boundaries, the cell
    centroids and the cell area
    """
    if 'z' not in spots.columns:
        # spots = spots.assign(z=np.zeros(spots.shape[0]))
        z = []
        n = len(coo)
        for i, row in spots.iterrows():
            r = divmod(i, n)[1]
            # r = uniform(0, n)
            z.append(r)
        spots = spots.assign(z=z)
        logger.info(' Added z-dimension to the spots dataframe')

    if isinstance(coo, list):
        for i in range(len(coo)):
            assert isinstance(coo[i], coo_matrix)

    if isinstance(coo, coo_matrix):
        assert coo.shape == 2
        coo = [coo]

    # 4. If a spots is placed withing the cell boundaries,
    # label it, otherwise the label = 0
    spots_df = label_spots(spots, coo)

    label_image_3d = np.stack([d.toarray().astype(np.uint32) for d in coo])
    logger.info(' Number of spots passed-in: %d' % spots.shape[0])
    logger.info(' Number of segmented cells: %d' % sum(np.unique(label_image_3d) > 0))
    # logger.info(' Segmentation array implies that image has width: %dpx and height: %dpx' % (coo.shape[1], coo.shape[0]))

    # 6. Get cell centroids and area
    properties = ['label', 'area', 'centroid', 'filled_image']
    props_df = pd.DataFrame(skmeas.regionprops_table(label_image_3d, properties=properties))
    num_slices = props_df.filled_image.apply(img_depth).values
    props_df = props_df.assign(mean_area_per_slice=props_df.area/num_slices)
    props_df = props_df.drop(['filled_image'], axis=1)
    props_df = props_df.rename(columns={"centroid-0": "z_cell", "centroid-1": "y_cell", "centroid-2": "x_cell"})

    # 7. Get the cell boundaries
    boundaries_list = []
    i = 0
    z_page = coo[i].toarray()
    b = extract_borders_dip(z_page.astype(np.uint32), 0, 0, [0])
    b = b.assign(z=i)
    b = b[['label', 'z', 'coords']]
    boundaries_list.append(b)
    cell_boundaries = pd.concat(boundaries_list)

    logger.info('Keeping boundaries for z=0 only')
    cell_boundaries = cell_boundaries[cell_boundaries.z == 0]  # Select the boundaries on z=0 only

    _cells = props_df.rename(columns={'x_cell': 'x', 'y_cell': 'y', 'z_cell': 'z'})
    _cell_boundaries = cell_boundaries.sort_values(['label', 'z'])
    _spots = spots_df[['x', 'y', 'z', 'label', 'Gene']].rename(columns={'Gene': 'target', 'x': 'x_global', 'y': 'y_global', 'z': 'z_global'})

    return _cells, _cell_boundaries, _spots


def img_depth(x):
    if len(x.shape) == 3:
        return x.shape[0]
    elif len(x.shape) == 2:
        return 1
    else:
        return np.nan

