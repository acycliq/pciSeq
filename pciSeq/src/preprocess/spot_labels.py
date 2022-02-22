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


# def inside_cell(label_image: np.array, idx: np.array) -> np.array:
#     """
#     :param label_image: An array of size height-by-width for the label image.
#     :param idx: An array of size 2-by-N of the pixels coordinates of spot idx[k], k=1...N
#     :return:
#     a = np.array([  [4,0,1],
#                     [2,0,0],
#                     [0,1,0]])
#
#     idx = np.array([[0,0],
#                     [2, 1],
#                     [1,2],
#                     [1,3]])
#
#     inside_cell(a, idx.T) = [4., 1., 0., nan]
#     which means that:
#             spot with coords [0,0] lies inside cell 4
#             spot with coords [2,0] lies inside cell 1
#             spot with coords [1,2] is a background spot
#             spot with coords [1,3] is outside the bounds and assigned to nan
#
#     """
#     assert isinstance(idx[0], np.ndarray), "Array 'idx' must be an array of arrays."
#     idx = idx.astype(np.int64)
#     out = np.array([])
#     dim = np.ones(idx.shape[0], dtype=int)
#     dim[:len(label_image.shape)] = label_image.shape
#
#     # output array
#     out = np.nan * np.ones(idx.shape[1], dtype=int)
#
#     # find the ones within bounds:
#     is_within = np.all(idx.T <= dim-1, axis=1)
#
#     # also keep only non-negative ones
#     is_positive = np.all(idx.T >= 0, axis=1)
#
#     # filter array
#     arr = idx[:, is_within & is_positive]
#     flat_idx = np.ravel_multi_index(arr, dims=dim, order='C')
#     out[is_within & is_positive] = label_image.ravel()[flat_idx]
#
#     # if the matrix a is a coo_matrix then the following should be
#     # equivalent (maybe better memory-wise since you do not have use
#     # a proper array (no need to do coo.toarray())
#     # out[is_within & is_positive] = a.tocsr(arr)
#     # print('in label_spot')
#
#     return out


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


def stage_data(spots: pd.DataFrame, coo: coo_matrix) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Reads the spots and the label image that are passed in and calculates which cell (if any) encircles any
    given spot within its boundaries. It also retrieves the coordinates of the cell boundaries, the cell
    centroids and the cell area
    """
    logger.info(' Number of spots passed-in: %d' % spots.shape[0])
    logger.info(' Number of segmented cells: %d' % len(set(coo.data)))
    logger.info(' Segmentation array implies that image has width: %dpx and height: %dpx' % (coo.shape[1], coo.shape[0]))
    mask_x = (spots.x >= 0) & (spots.x <= coo.shape[1])
    mask_y = (spots.y >= 0) & (spots.y <= coo.shape[0])
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
    inc = inside_cell(coo.tocsr(), spots)
    spots = spots.assign(label=inc)

    # 2. Get cell centroids and area
    props = skmeas.regionprops(coo.toarray().astype(np.int32))
    props_df = pd.DataFrame(data=[(d.label, d.area, d.centroid[1], d.centroid[0]) for d in props],
                      columns=['label', 'area', 'x_cell', 'y_cell'])

    # 3. Get the cell boundaries
    cell_boundaries = extract_borders_dip(coo.toarray().astype(np.uint32), 0, 0, [0])

    assert props_df.shape[0] == cell_boundaries.shape[0] == coo.data.max()
    assert set(spots.label[spots.label > 0]) <= set(props_df.label)

    cells = props_df.merge(cell_boundaries)
    cells.sort_values(by=['label', 'x_cell', 'y_cell'])
    assert cells.shape[0] == cell_boundaries.shape[0] == props_df.shape[0]

    # join spots and cells on the cell label so you can get the x,y coords of the cell for any given spot
    spots = spots.merge(cells, how='left', on=['label'])

    _cells = cells[['label', 'area', 'x_cell', 'y_cell']].rename(columns={'x_cell': 'x', 'y_cell': 'y'})
    _cell_boundaries = cells[['label', 'coords']]
    _spots = spots[['x', 'y', 'label', 'Gene', 'x_cell', 'y_cell']].rename(columns={'Gene': 'target', 'x': 'x_global', 'y': 'y_global'})

    return _cells, _cell_boundaries, _spots
