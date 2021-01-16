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
from scipy.sparse import coo_matrix, save_npz, load_npz
from pciSeq.src.preprocess.cell_borders import extract_borders_par, extract_borders_dip
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()


def inside_cell(a, idx):
    '''
    Given an array a (image_array) and
    :param a: An array of size numPixelsY-by-numPixelsX specifying that element (i,j) belongs to
                cell a[i,j]. Note that array a is 1-based, ie if pixel (i,j) is outside a cell then
                a[i,j] = 0.
    :param idx: An array of size 2-by-N of the pixels coordinates of spot idx[k], k=1...N
    :return:
    a = np.array([  [4,0,1],
                    [2,0,0],
                    [0,1,0]])

    idx = np.array([[0,0],
                    [2, 1],
                    [1,2],
                    [1,3]])

    IndexArrayNan(a, idx.T) = [4., 1., 0., nan]
    which means that:
            spot with coords [0,0] belongs to cell 4
            spot with coords [2,0] belongs to cell 1
            spot with coords [1,2] belongs to 0 (ie no assigned cell)
            spot with coords [1,3] is outside the bounds and assigned to nan

    '''
    assert isinstance(idx[0], np.ndarray), "Array 'idx' must be an array of arrays."
    idx = idx.astype(np.int64)
    out = np.array([])
    dim = np.ones(idx.shape[0], dtype=int)
    dim[:len(a.shape)] = a.shape

    # output array
    out = np.nan * np.ones(idx.shape[1], dtype=int)

    # find the ones within bounds:
    is_within = np.all(idx.T <= dim-1, axis=1)

    # also keep only non-negative ones
    is_positive = np.all(idx.T >= 0, axis=1)

    # filter array`
    arr = idx[:, is_within & is_positive]
    flat_idx = np.ravel_multi_index(arr, dims=dim, order='C')
    out[is_within & is_positive] = a.ravel()[flat_idx]

    # if the matrix a is a coo_matrix then the following should be
    # equivalent (maybe better memory-wise since you do not have use
    # a proper array (no need to do coo.toarray())
    # out[is_within & is_positive] = a.tocsr(arr)
    # print('in label_spot')

    return out


def remap_labels(coo):
    """
    Used for debugging only. It resuffles the label_image
    """
    coo_max = coo.data.max()
    _keys = 1 + np.arange(coo_max)
    _vals = _keys.copy()
    np.random.shuffle(_vals)
    d = dict(zip(_keys, _vals))
    new_data = np.array([d[x] for x in coo.data]).astype(np.uint64)
    out = coo_matrix((new_data, (coo.row, coo.col)), shape=coo.shape)
    return out


def stage_data(spots, coo):
    logger.info('Number of spots passed in: %d' % spots.shape[0])
    logger.info('Number of segmented cells %d' % len(set(coo.data)))
    logger.info('Segmentation array implies that image has width: %dpx and height: %dpx' % (coo.shape[1], coo.shape[0]))
    # spots = pd.read_csv(cfg['spots'])
    # coo = load_npz(cfg['label_image'])
    # spots_df = spots_df[['Gene', 'xc', 'yc']].rename(columns={'xc': 'x', 'yc': 'y'})
    # spots_df.x = spots_df.x - 6150
    # spots_df.y = spots_df.y - 12987
    mask_x = (spots.x >= 0) & (spots.x <= coo.shape[1])
    mask_y = (spots.y >= 0) & (spots.y <= coo.shape[0])
    spots = spots[mask_x & mask_y]

    # resuffle
    # spots = spots.sample(frac=1).reset_index(drop=True)

    # _point = [5471-14, 110]
    # logger.info('label at (y, x): (%d, %d) is %d' % (_point[0], _point[1], coo.toarray()[_point[0], _point[1]]))

    # coo = remap_labels(coo)
    # logger.info('remapped label at (y, x): (%d, %d) is %d' % (_point[0], _point[1], coo.toarray()[_point[0], _point[1]]))

    yx_coords = spots[['y', 'x']].values.T
    inc = inside_cell(coo.toarray(), yx_coords)
    spots = spots.assign(label=inc)

    props = skmeas.regionprops(coo.toarray().astype(np.int32))
    props_df = pd.DataFrame(data=[(d.label, d.area, d.centroid[1], d.centroid[0]) for d in props],
                      columns=['label', 'area', 'x_cell', 'y_cell'])

    cell_boundaries = extract_borders_dip(coo.toarray(), 0, 0, [0])

    assert props_df.shape[0] == cell_boundaries.shape[0] == coo.data.max()
    assert set(spots.label[spots.label > 0]) <= set(props_df.label)

    cells = props_df.merge(cell_boundaries)
    cells.sort_values(by=['label', 'x_cell', 'y_cell'])
    assert cells.shape[0] == cell_boundaries.shape[0] == props_df.shape[0]

    spots = spots.merge(cells, how='left', on=['label'])

    # dirpath = cfg['temp']
    # _cells, _cell_boundaries, _spots = writer(cells, spots, dirpath)
    # return _cells, _cell_boundaries, _spots
    _cells = cells[['label', 'area', 'x_cell', 'y_cell']].rename(columns={'x_cell': 'x', 'y_cell': 'y'})
    _cell_boundaries = cells[['label', 'coords']]
    _spots = spots[['x', 'y', 'label', 'Gene', 'x_cell', 'y_cell']].rename(columns={'Gene': 'target', 'x': 'x_global', 'y': 'y_global'})
    return _cells, _cell_boundaries, _spots


def writer(cells, spots, dirpath):
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    _cells, _cell_boundaries = writer_cells(cells, dirpath)
    _spots = writer_spots(spots, dirpath)
    return _cells, _cell_boundaries, _spots


def writer_cells(cell_props, dirpath):
    '''
    save the data to the flatfile
    :return:
    '''

    cell_props = cell_props.sort_values(by=['label', 'x_cell', 'y_cell'])

    # 1. save the cell props
    cell_props = cell_props.rename({'x_cell': 'x', 'y_cell': 'y'}, axis=1)
    # cell_props['x'] = cell_props.x.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
    # cell_props['y'] = cell_props.y.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
    cell_props['cell_id'] = cell_props.cell_id.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
    cell_props['label'] = cell_props.label.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

    cells_headers = ['cell_id', 'label', 'area', 'x', 'y']
    cell_props[cells_headers].to_csv(os.path.join(dirpath, '_cells.csv'), index=False)

    # 2. save the cell coords
    coords_headers = ['cell_id', 'label', 'coords']
    cell_props[coords_headers].to_json(os.path.join(dirpath, '_cellCoords.json'), orient='records')
    return cell_props[cells_headers], cell_props[coords_headers]


def writer_spots(spots_df, dirpath):

    # spots_df = spots_df.sort_values(by=['label', 'x', 'y'])

    # 3. save the spots
    spots_df['target'] = spots_df.Gene
    spots_df['x_global'] = spots_df.x
    spots_df['y_global'] = spots_df.y
    # spots_df['x_cell'] = spots_df.x_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
    # spots_df['y_cell'] = spots_df.y_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

    spots_headers = ['x_global', 'y_global', 'label', 'target', 'x_cell', 'y_cell']
    spots_df[spots_headers].to_csv(os.path.join(dirpath, '_spots.csv'), index=False)
    print('Total number of collected spots: %d' % spots_df.shape[0])

    return spots_df[spots_headers]

### FUNCIONS BELOW ARE DEPRECATED ###
def coofy(spots, label_image):
# FUNCTION DEPRECATED
    x = spots.x.values.astype(int)
    y = spots.y.values.astype(int)

    coo_shape = label_image.shape
    idx = spots.index.values + 1  # avoid having idx = 0
    coo = coo_matrix((idx, (y, x)), shape=coo_shape)
    return coo


def _spot_parent_label(spots, sp_label_image):
# DEPRECATED. Replaced by inside_cell()

    if sp_label_image.nnz == 0:
        # The tile empty, all spots are on the background
        spots['label'] = 0
    else:
        # 1. get the label_image
        label_image = sp_label_image.toarray()

        # 2. unscaled and local coordinates for the spots (as a sparse array)
        coo = coofy(spots.copy(), label_image)

        coo_arr = coo.toarray()
        label_coords = coo_matrix(label_image * coo_arr.astype(bool))
        spot_id_coords = coo_matrix(label_image.astype(bool) * coo_arr)

        df = pd.DataFrame({'x': label_coords.col,
                           'y': label_coords.row,
                           'label': label_coords.data},
                          index=spot_id_coords.data).astype(int)

        if np.any(df.index.duplicated()):
            print('Found %d duplicated. Investigate!' % df.index.duplicated().sum())
            ## that actually means that the same spot (ie read/dot) exists at two different locations at the same time

        df = df[~df.index.duplicated()]
        # spots['label'] = df.label
        spots = spots.merge(df, how='left', on=['x', 'y'])

        # if nan it means the spot is on the background. Hence set the label = 0
        spots['label'] = spots.label.fillna(0).astype(int)
    return spots

