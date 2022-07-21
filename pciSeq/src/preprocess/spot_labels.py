"""
Functions to prepare the data for pciSeq. The label image and spots are parsed and if a spot
lies within the cell boundaries then the corresponding cell id is recorded.
Cell centroids and cell areas are also calculated.
"""

import numpy as np
import pandas as pd
from typing import Tuple
from scipy.sparse import coo_matrix, csr_matrix
from pciSeq.src.preprocess.cell_borders import extract_borders_dip
from pciSeq.src.preprocess.regionprops import regionprops
from pciSeq.src.cell_call.log_config import logger


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

    labels = np.concatenate([np.unique(d.data) for d in label_image])
    if labels.max() != len(set(labels)):
        logger.info(' Labels in segmentation array are not a sequence of integers without gaps between them. Relabelling...')
        labels = np.append(0, labels) # append 0 (the background label)
        _, idx = np.unique(labels, return_inverse=True)
        assert idx[0] == 0  # make sure the background is at position 0
        dic = dict(zip(labels, idx.astype(np.uint32)))
        assert dic[0] == 0, "The background label (ie 0) should be the smallest of the labels"
        out = []
        for plane in label_image:
            row = plane.row
            col = plane.col
            data = [dic[d] for d in plane.data]
            shape = plane.shape
            c = coo_matrix((data, (row, col)), shape=shape)
            out.append(c)
    else:
        out = label_image
    return out


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


def weighted_average(df,data_col,weight_col,by_col):
    df['_data_times_weight'] = df[data_col]*df[weight_col]
    df['_weight_where_notnull'] = df[weight_col]*pd.notnull(df[data_col])
    g = df.groupby(by_col)
    result = g['_data_times_weight'].sum() / g['_weight_where_notnull'].sum()
    del df['_data_times_weight'], df['_weight_where_notnull']
    return result


def stage_data(spots: pd.DataFrame, coo: coo_matrix, cfg) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Reads the spots and the label image that are passed in and calculates which cell (if any) encircles any
    given spot within its boundaries. It also retrieves the coordinates of the cell boundaries, the cell
    centroids and the cell area
    """
    if cfg['is_3D'] or cfg['relax_segmentation']:
        z_min = cfg['3D:from_plane_num']
        z_max = cfg['3D:to_plane_num']
        spots = spots.assign(z=spots.z_stack * cfg['3D:anisotropy'])
        if cfg['is_3D']:
            spots, coo = truncate_data(spots, coo, z_min, z_max)
            coo = remove_cells(coo)
    else:
        spots = spots.assign(z=spots.z_stack)

    coo = reorder_labels(coo)

    img_shape = set([d.shape for d in coo])
    assert len(img_shape) == 1, 'pages do not have the same shape'
    img_shape = img_shape.pop()
    logger.info(' Number of spots passed-in: %d' % spots.shape[0])
    logger.info(' Number of segmented cells: %d' % max([d.data.max() for d in coo if len(d.data) > 0]))
    logger.info(' Segmentation array implies that image has width: %dpx and height: %dpx' % (list(img_shape)[1], list(img_shape)[0]))
    mask_x = (spots.x >= 0) & (spots.x <= img_shape[1])
    mask_y = (spots.y >= 0) & (spots.y <= img_shape[0])
    spots = spots[mask_x & mask_y]

    # 1. Find which cell the spots lie within
    spots = spots.assign(label=np.zeros(spots.shape[0]))
    for z in np.unique(spots.z_stack):
        spots_z = spots[spots.z_stack == z]
        inc = inside_cell(coo[int(z)].tocsr(), spots_z)
        spots.loc[spots.z_stack == z, 'label'] = inc

    # 2. Get cell centroids and area
    properties = ("label", "area", "centroid")
    props = regionprops(coo, properties=properties).compute()

    centroid_0 = weighted_average(props, 'dim-0', 'area', 'label')
    y_cell = weighted_average(props, 'centroid-0', 'area', 'label')
    x_cell = weighted_average(props, 'centroid-1', 'area', 'label')
    mean_area_per_slice = weighted_average(props, 'area', 'area', 'label')

    assert np.all(centroid_0.index == y_cell.index)
    assert np.all(centroid_0.index == x_cell.index)
    assert np.all(centroid_0.index == mean_area_per_slice.index)
    assert centroid_0.index.max() == centroid_0.shape[0]
    props_df = pd.DataFrame({
        'label': centroid_0.index,
        'centroid-0':  centroid_0.values,
        'y_cell': y_cell.values,
        'x_cell': x_cell.values,
        'z_cell': centroid_0.values * cfg['3D:anisotropy'],
        'area': mean_area_per_slice.values, # mean area per slice
    })

    # 3. Get the cell boundaries, Only needed for the 2D case
    if not cfg['is_3D'] or cfg['relax_segmentation']:
        _cell_boundaries = extract_borders_dip(coo[0].toarray().astype(np.uint32), 0, 0, [0])
        # If you relax the segmentation constraint then do the ellipsoid borders as well
        if cfg['relax_segmentation']:
           pass
    else:
        _cell_boundaries = None

    _cells = props_df.rename(columns={'x_cell': 'x', 'y_cell': 'y', 'z_cell': 'z'})
    _spots = spots[['x', 'y', 'z', 'label', 'Gene']]
    _cell_boundaries = _cell_boundaries.rename(columns={'label': 'cell_id'})

    return _cells, _cell_boundaries, _spots

