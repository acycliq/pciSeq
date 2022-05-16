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


def reorder_labels(coo):
    """
    rearranges the labels so that they are a sequence of integers
    """
    label_image = np.stack([d.toarray().astype(np.uint16) for d in coo])
    _, idx = np.unique(label_image.flatten(), return_inverse=True)
    label_image = idx.reshape(label_image.shape)
    ## Need to check what is happening when you change to no.uint16. For example coo[-1].data.sum() is not the same
    ## when you change to np.uint16
    return [coo_matrix(d.astype(np.uint64)) for d in label_image]


def remove_cells(coo_list):
    """
    removes cells that exist on just one single frame of the 
    z-stack
    """
    labels_per_frame = [np.unique(d.data) for d in coo_list]
    label_counts = np.bincount([d for labels in labels_per_frame for d in labels])
    single_page_labels = [d[0] for d in enumerate(label_counts) if d[1] == 1]
    removed_cells = []
    _frames = []
    for i, coo in enumerate(coo_list):
        s = set(coo.data).intersection(set(single_page_labels))
        for d in s:
            coo.data[coo.data == d] = 0
            removed_cells.append(d)
            # logger.info('Removed cell:%d from frame: %d' % (d, i))
        _frames.append(i)
    len_c = len(set(removed_cells))
    len_f = len(set(_frames))
    logger.info(' Found %d cells that exist on just one single frame. Those cells have been removed from %i frames.' % (len_c, len_f))
    return coo_list


def truncate_data(spots, label_image, z_min, z_max):
    if z_min is not None and z_max is not None:
        logger.info(' Truncating masks and spots. Keeping those between frame %d and %d only' % (z_min, z_max))
        spots = truncate_spots(spots, z_min, z_max)
        label_image = truncate_zstack(label_image, z_min, z_max)
    return spots, label_image


def truncate_zstack(masks, z_min, z_max):
    return masks[z_min: z_max+1]


def truncate_spots(spots, zmin, zmax):
    spots_min = spots[(spots.z_stack <= zmax) & (spots.z_stack >= zmin)]
    spots_min = spots_min.assign(z_stack=spots_min.z_stack - zmin)
    # out = spots_min.z - zmin
    return spots_min


def stage_data(spots: pd.DataFrame, coo: coo_matrix, cfg) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Reads the spots and the label image that are passed in and calculates which cell (if any) encircles any
    given spot within its boundaries. It also retrieves the coordinates of the cell boundaries, the cell
    centroids and the cell area
    """
    ppm = cfg['ppm']
    z_min = cfg['z_stack_min']
    z_max = cfg['z_stack_max']
    # z_min = 18
    # z_max = 43
    spots, coo = truncate_data(spots, coo, z_min, z_max)
    coo = remove_cells(coo)
    coo = reorder_labels(coo)

    label_image = np.stack([d.toarray().astype(np.uint16) for d in coo])
    logger.info(' Number of spots passed-in: %d' % spots.shape[0])
    logger.info(' Number of segmented cells: %d' % sum(np.unique(label_image) > 0) )
    logger.info(' Segmentation array implies that image has width: %dpx and height: %dpx' % (label_image.shape[2], label_image.shape[1]))
    mask_x = (spots.x >= 0) & (spots.x <= label_image.shape[2])
    mask_y = (spots.y >= 0) & (spots.y <= label_image.shape[1])
    spots = spots[mask_x & mask_y]

    # Debugging code!
    # resuffle
    # spots = spots.sample(frac=1).reset_index(drop=True)

    # _point = [5471-14, 110]
    # logger.info('label at (y, x): (%d, %d) is %d' % (_point[0], _point[1], coo.toarray()[_point[0], _point[1]]))

    # coo = remap_labels(coo)
    # logger.info('remapped label at (y, x): (%d, %d) is %d' % (_point[0], _point[1], coo.toarray()[_point[0], _point[1]]))

    # 1. Find which cell the spots lie within
    # yx_coords = spots[['y', 'x']].values.T
    spots = spots.assign(label=np.zeros(spots.shape[0]))
    for z in np.unique(spots.z_stack):
        spots_z = spots[spots.z_stack == z]
        inc = inside_cell(coo[int(z)].tocsr(), spots_z)
        spots.loc[spots.z_stack == z, 'label'] = inc
        # spots = spots.assign(label=inc)

    # 2. Get cell centroids and area
    properties = ['label', 'area', 'centroid', 'filled_image']
    props_df = pd.DataFrame(skmeas.regionprops_table(label_image, properties=properties))
    num_slices = props_df.filled_image.apply(img_depth).values
    props_df = props_df.assign(mean_area_per_slice=props_df.area/num_slices)
    props_df = props_df.assign(z_cell=props_df['centroid-0'] * ppm)
    props_df = props_df.drop(['filled_image'], axis=1)
    props_df = props_df.rename(columns={"centroid-1": "y_cell", "centroid-2": "x_cell"})

    # 7. Get the cell boundaries
    boundaries_list = []
    # for i, d in enumerate(coo):
    i = 0
    d = coo[0]
    z_page = d.toarray()
    b = extract_borders_dip(z_page.astype(np.uint32), 0, 0, [0])
    b = b.assign(z=i)
    b = b[['label', 'z', 'coords']]
    boundaries_list.append(b)
    cell_boundaries = pd.concat(boundaries_list)

    # logger.info('Keeping boundaries for z=0 only')
    cell_boundaries = cell_boundaries[cell_boundaries.z == 0]  # Select the boundaries on z=0 only

    _cells = props_df.rename(columns={'x_cell': 'x', 'y_cell': 'y', 'z_cell': 'z'})
    _cell_boundaries = cell_boundaries.sort_values(['label', 'z'])
    _spots = spots[['x', 'y', 'z', 'label', 'Gene']].rename(columns={'Gene': 'target', 'x': 'x_global', 'y': 'y_global', 'z': 'z_global'})

    return _cells, _cell_boundaries, _spots


def img_depth(x):
    if len(x.shape) == 3:
        return x.shape[0]
    elif len(x.shape) == 2:
        return 1
    else:
        return np.nan
