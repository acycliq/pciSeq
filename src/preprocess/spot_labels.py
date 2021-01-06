import os
import shutil
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix, save_npz, load_npz
import skimage.measure as skmeas
from src.preprocess.cell_borders import extract_borders_par, extract_borders_dip
import logging


dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()

def coofy(spots, label_image):
    x = spots.x.values.astype(int)
    y = spots.y.values.astype(int)

    coo_shape = label_image.shape
    idx = spots.index.values + 1  # avoid having idx = 0
    coo = coo_matrix((idx, (y, x)), shape=coo_shape)
    return coo


def _spot_parent_label(spots, sp_label_image):
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


def writer(cells, spots, dirpath):
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    writer_cells(cells, dirpath)
    writer_spots(spots, dirpath)


def writer_cells(cell_props, dirpath):
    '''
    save the data to the flatfile
    :return:
    '''

    cell_props = cell_props.sort_values(by=['label', 'x_cell', 'y_cell'])

    # 1. save the cell props
    cell_props = cell_props.rename({'x_cell': 'x', 'y_cell': 'y'}, axis=1)
    cell_props['x'] = cell_props.x.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
    cell_props['y'] = cell_props.y.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
    cell_props['cell_id'] = cell_props.cell_id.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
    cell_props['label'] = cell_props.label.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

    cells_headers = ['cell_id', 'label', 'area', 'x', 'y']
    cell_props[cells_headers].to_csv(os.path.join(dirpath, '_cells.csv'), index=False)

    # 2. save the cell coords
    coords_headers = ['cell_id', 'label', 'coords']
    cell_props[coords_headers].to_json(os.path.join(dirpath, '_cellCoords.json'), orient='records')
    return cell_props[cells_headers], cell_props[coords_headers]


def writer_spots(spots_df, dirpath):

    spots_df = spots_df.sort_values(by=['label', 'x', 'y'])

    # 3. save the spots
    spots_df['target'] = spots_df.Gene
    spots_df['x_global'] = spots_df.x
    spots_df['y_global'] = spots_df.y
    spots_df['x_cell'] = spots_df.x_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
    spots_df['y_cell'] = spots_df.y_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

    spots_headers = ['x_global', 'y_global', 'label', 'target', 'x_cell', 'y_cell']
    spots_df[spots_headers].to_csv(os.path.join(dirpath, '_spots.csv'), index=False)
    print('Total number of collected spots: %d' % spots_df.shape[0])

    return spots_df[spots_headers]


def remap_labels(coo):
    coo_max = coo.data.max()
    _keys = 1 + np.arange(coo_max)
    _vals = _keys.copy()
    np.random.shuffle(_vals)
    d = dict(zip(_keys, _vals))
    new_data = np.array([d[x] for x in coo.data]).astype(np.uint64)
    out = coo_matrix((new_data, (coo.row, coo.col)), shape=coo.shape)
    return out


def stage_data(cfg):
    spots = pd.read_csv(cfg['spots'])
    coo = load_npz(cfg['label_image'])
    # spots_df = spots_df[['Gene', 'xc', 'yc']].rename(columns={'xc': 'x', 'yc': 'y'})
    # spots_df.x = spots_df.x - 6150
    # spots_df.y = spots_df.y - 12987
    mask_x = (spots.x >= 0) & (spots.x <= coo.shape[1])
    mask_y = (spots.y >= 0) & (spots.y <= coo.shape[0])
    spots = spots[mask_x & mask_y]

    # resuffle
    # spots = spots.sample(frac=1).reset_index(drop=True)
    _point = [14, 110]
    logger.info('label at (y, x): (%d, %d) is %d' % (_point[0], _point[1], coo.toarray()[_point[0], _point[1]]))
    coo = remap_labels(coo)
    logger.info('remapped label at (y, x): (%d, %d) is %d' % (_point[0], _point[1], coo.toarray()[_point[0], _point[1]]))



    spot_label = _spot_parent_label(spots, coo)

    props = skmeas.regionprops(coo.toarray().astype(np.int32))
    props_df = pd.DataFrame(data=[(d.label, d.area, d.centroid[1], d.centroid[0]) for d in props],
                      columns=['label', 'area', 'x_cell', 'y_cell'])

    cell_boundaries = extract_borders_dip(coo.toarray(), 0, 0, [0])

    assert props_df.shape[0] == cell_boundaries.shape[0] == coo.data.max()
    assert set(spot_label.label[spot_label.label > 0]) <= set(props_df.label)

    cells = props_df.merge(cell_boundaries)
    cells['cell_id'] = cells.label - 1
    cells.sort_values(by=['label', 'x_cell', 'y_cell'])
    assert cells.shape[0] == cell_boundaries.shape[0] == props_df.shape[0]

    spots = spot_label.merge(cells, how='left', on=['label'])
    spots.sort_values(by=['label', 'x', 'y'])
    assert spots.shape[0] == spot_label.shape[0]

    dirpath = cfg['temp']
    writer(cells, spots, dirpath)


if __name__ == "__main__":
    # spots = pd.read_csv('data/mouse/ca1/iss/spots.csv')
    # spots_df = pd.read_csv('D:\\Home\\Dimitris\\OneDrive - University College London\\dev\\Python\\cell_call\\demo_data\\spots_sanity_check.csv')

    spots_df = pd.read_csv('../../data/mouse/ca1/iss/spots.csv')
    label_image = load_npz('../../data/mouse/ca1/segmentation/label_image.coo.npz')
    spot_labels(spots_df, label_image)

    print('Done')
