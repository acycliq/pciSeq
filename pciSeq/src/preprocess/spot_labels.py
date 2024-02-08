"""
Functions to prepare the data for pciSeq. The label image and spots are parsed and if a spot
lies within the cell boundaries then the corresponding cell id is recorded.
Cell centroids and cell areas are also calculated.
"""

import numpy as np
import pandas as pd
import skimage.measure as skmeas
from typing import Tuple
from multiprocessing.dummy import Pool as ThreadPool
from scipy.sparse import coo_matrix, csr_matrix
from pciSeq.src.preprocess.cell_borders import extract_borders_dip, extract_borders
from pciSeq.src.core.utils import get_img_shape, adjust_for_anisotropy
import logging

spot_labels_logger = logging.getLogger(__name__)


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


def get_unique_labels(masks):
    """
    same as:
        [np.unique(d.data) for d in label_image]
    but faster
    """
    pool = ThreadPool()  # ThreadPool??!! You sure? maybe process pool
    out = pool.map(lambda mask: np.unique(mask.data), masks)
    pool.close()
    pool.join()
    return out


def reorder_labels(label_image):
    """
    rearranges the labels so that they are a sequence of integers

    Parameters
    ----------
    label_image: A 2D, 3D array or a list of coo_matrices. If a 3D array, the first dimension is
                the Z-axis, the axis to loop over to get the planes of the image

    Returns
    ------
    out:    a list of coo_matrices. Each element of the list is a sparse matrix representing a plane
            of the segmented image. The labels have been relabelled (if needed and across the full stack) so that
            they form a sequence of integers
    """
    try:
        assert isinstance(label_image, np.ndarray)
        assert len(label_image.shape) in {2, 3}
        label_image = [coo_matrix(d) for d in label_image]
    except AssertionError:
        assert np.all([isinstance(d, coo_matrix) for d in label_image]) and isinstance(label_image, list)
    else:
        raise TypeError('input label_image should be an 2D, 3D array or a list of coo matrices')

    # labels = np.concatenate([np.unique(d.data) for d in label_image])
    labels = np.concatenate(get_unique_labels(label_image))
    if labels.max() != len(set(labels)):
        spot_labels_logger.info(' Labels in segmentation array are not a sequence of integers without gaps between them. Relabelling...')
        labels = np.append(0, labels)  # append 0 (the background label). Makes sure the background is at position 0
        _, idx = np.unique(labels, return_inverse=True)
        label_dict = dict(zip(labels, idx.astype(np.uint32))) # maps the cellpose ids to the new ids
        assert label_dict[0] == 0, "The background label (ie 0) should be the smallest of the labels"
        out = []
        for plane in label_image:
            row = plane.row
            col = plane.col
            data = [label_dict[d] for d in plane.data]
            shape = plane.shape
            c = coo_matrix((data, (row, col)), shape=shape)
            out.append(c)

        label_map = pd.DataFrame(label_dict.items(), columns=['old_label', 'new_label'])
    else:
        out = label_image
        label_map = None
    return out, label_map


def truncate_label_image(coo, cfg):
    arr = np.arange(len(coo))
    coo = [coo[d] for d in arr if d not in cfg['exclude_planes']]
    return coo


def truncate_spots(spots, cfg):
    int_z = np.floor(spots.z_plane)
    mask = [True if d not in cfg['exclude_planes'] else False for d in int_z]
    spots = spots[mask]

    # find the first plane we should be keeping. Assumes that exclude_list is
    # in increasing order. Made sure that is true in the validation step
    diff = np.diff(cfg['exclude_planes']) - 1
    if np.all(diff == 0):
        # we remove planes from the bottom of the z-stack only,
        # hence the first plane that we keep is the first one
        # that is not in the exclude list
        min_plane = max(cfg['exclude_planes']) + 1
    else:
        # otherwise we remove plane from the bottom and the top of
        # the z-stack, possibly (but unlikely) from the middle too.
        # Start scanning the z-stack from the bottom, find where the
        # discontinuity happens and get the first plane we keep
        iLeft = list(diff > 0).index(True)
        min_plane = cfg['exclude_planes'][iLeft] + 1

    spots.loc[:, 'z_plane'] = spots.z_plane - min_plane
    return spots, min_plane


def remove_cells(coo_list, cfg):
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
            # _frames.append(i + cfg['from_plane_num'])
            _frames.append(i)
            # logger.info('Removed cell:%d from frame: %d' % (d, i))
    len_c = len(set(removed_cells))
    len_f = len(set(_frames))

    if len(removed_cells) > 0:
        spot_labels_logger.info('Found %d cells that exist on just one single plane. Those cells have been removed '
                                'from %i planes.' % (len_c, len_f))

    removed_df = pd.DataFrame({
        'removed_cell_label': removed_cells,
        'frame_num': _frames,
        'comment': 'labels are the original labels as they appear in the segmentation masks that were passed-in to '
                   'pciSeq'
    })
    return coo_list, removed_df


def remove_oob(spots, img_shape):
    """
    removes out of bounds spots (if any...)
    """
    mask_x = (spots.x >= 0) & (spots.x <= img_shape[1])
    mask_y = (spots.y >= 0) & (spots.y <= img_shape[0])
    return spots[mask_x & mask_y]


def stage_data(spots: pd.DataFrame, coo: coo_matrix, cfg: dict) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Reads the spots and the label image that are passed in and calculates which cell (if any) contains any
    given spot within its boundaries. It also retrieves the coordinates of the cell boundaries, the cell
    centroids and the cell area
    """

    # drop planes in the exclude planes list
    if cfg['is3D'] and cfg['exclude_planes'] is not None:
        coo = truncate_label_image(coo, cfg)
        spots, min_plane = truncate_spots(spots, cfg)
        coo, removed = remove_cells(coo, cfg)
        if removed.shape[0] > 0:
            # add the min_plane back so that the planes refer to the original 3d image before
            # the the removal of the planes described in the exclude_planes list
            removed.frame_num = removed.frame_num + min_plane

    coo, label_map = reorder_labels(coo)
    [h, w] = get_img_shape(coo)
    spots = remove_oob(spots.copy(), [h, w])
    spots_xyz = adjust_for_anisotropy(spots, cfg['voxel_size'])


    label_map = None
    if coo.data.max() != len(set(coo.data)):
        spot_labels_logger.info('The labels in the label image do not seem to be a sequence of successive integers. Relabelling the label image.')
        coo, label_map = reorder_labels(coo)

    spot_labels_logger.info('Number of spots passed-in: %d' % spots.shape[0])
    spot_labels_logger.info('Number of segmented cells: %d' % len(set(coo.data)))
    spot_labels_logger.info('Segmentation array implies that image has width: %dpx and height: %dpx' % (coo.shape[1], coo.shape[0]))
    mask_x = (spots.x >= 0) & (spots.x <= coo.shape[1])
    mask_y = (spots.y >= 0) & (spots.y <= coo.shape[0])
    spots = spots[mask_x & mask_y]

    # 1. Find which cell the spots lie within
    inc = inside_cell(coo.tocsr(), spots)
    spots = spots.assign(label=inc)

    # 2. Get cell centroids and area
    props = skmeas.regionprops_table(coo.toarray().astype(np.int32), properties=['label', 'area', 'centroid'])
    props_df = pd.DataFrame(props).rename(columns={'centroid-0': 'y_cell', 'centroid-1': 'x_cell'})

    # if there is a label map, attach it to the cell props.
    if label_map is not None:
        props_df = pd.merge(props_df, label_map, left_on='label', right_on='new_label', how='left')
        props_df = props_df.drop(['new_label'], axis=1)

    # 3. Get the cell boundaries
    # cell_boundaries = extract_borders_dip(coo.toarray().astype(np.uint32))
    cell_boundaries = extract_borders(coo.toarray().astype(np.uint32))
    assert props_df.shape[0] == cell_boundaries.shape[0] == np.unique(coo.data).shape[0]
    assert set(spots.label[spots.label > 0]) <= set(props_df.label)

    cells = props_df.merge(cell_boundaries)
    cells.sort_values(by=['label', 'x_cell', 'y_cell'])
    assert cells.shape[0] == cell_boundaries.shape[0] == props_df.shape[0]

    # join spots and cells on the cell label so you can get the x,y coords of the cell for any given spot
    spots = spots.merge(cells, how='left', on=['label'])

    _cells = cells.drop(columns=['coords'])
    _cells = _cells.rename(columns={'x_cell': 'x0', 'y_cell': 'y0'})
    _cell_boundaries = cells[['label', 'coords']].rename(columns={'label': 'cell_id'})
    _spots = spots[['x', 'y', 'label', 'Gene', 'x_cell', 'y_cell']].rename(columns={'Gene': 'target', 'x': 'x_global', 'y': 'y_global'})

    return _cells, _cell_boundaries, _spots

