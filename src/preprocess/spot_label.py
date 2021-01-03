"""
Functions to label a spot, ie to give each spot a label which reflects the cell whose soma covers the coords of the spots.
If a spot is plotted against the background, then label = 0
"""
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

logger = logging.getLogger()


def spot_label(stage):
    res_list = []
    for i, tile in enumerate(stage.tiles):
        assert i == tile['tile_id'], 'The list is miss-aligned'
        # Assign the parent cell and its coordinates
        res = _assign_spot_parent_helper(stage, tile)
        # tile['spots'] = res
        res['tile_id'] = tile['tile_id']
        res_list.append(res)
    return pd.concat(res_list).astype({"label": int})


def _assign_spot_parent_helper(stage, tile):
    '''
    assign the parent cell
    :return:
    '''

    # first find the label of the parent cell
    spots_df = tile['spots']
    label_image = tile['label_image']
    spots_temp = _spot_parent_label(tile, spots_df, label_image)

    # find now the coordinates of the parent cell
    out = spot_parent_coords(stage, spots_temp.copy())
    return out


def _spot_parent_label(tile, spots, sp_label_image):
    # spots = tile['spots']
    if sp_label_image.nnz == 0:
        # The tile empty, all spots are on the background
        spots['label'] = 0
    else:
        # 1. get the label_image
        label_image = sp_label_image.toarray()

        # 2. unscaled and local coordinates for the spots (as a sparse array)
        coo = coofy(spots.copy(), tile)

        coo_arr = coo.toarray()
        label_coords = coo_matrix(label_image * coo_arr.astype(bool))
        spot_id_coords = coo_matrix(label_image.astype(bool) * coo_arr)

        df = pd.DataFrame({'x': label_coords.col,
                           'y': label_coords.row,
                           'label': label_coords.data},
                          index=spot_id_coords.data).astype(int)

        if np.any(df.index.duplicated()):
            logger.warning('Found %d duplicated. Investigate!' % df.index.duplicated().sum())
            ## that actually means that the same spot (ie read/dot) exists at two different locations at the same time

        df = df[~df.index.duplicated()]
        # spots['label'] = df.label
        spots = spots.merge(df, how='left', on=['x', 'y'])

        # if nan it means the spot is on the background. Hence set the label = 0
        spots['label'] = spots.label.fillna(0).astype(int)
    return spots


def spot_parent_coords(stage, spots):
    cells = stage.cell_props[['label', 'x', 'y']].rename(columns={'x': 'x_cell', 'y': 'y_cell'})
    out = spots.merge(cells, how='left', on=['label'])
    out['label'] = out.label.astype('Int64')
    return out


def localise_coords(spots, tile):
    '''
    convert spots coords to local coords
    (lacal means the origin is the top left corner of the tile not of the full image)
    :param spots:
    :return:
    '''
    # spots = tile['spots']
    df = spots[['x', 'y']]
    x0 = tile['tile_offset_x']
    y0 = tile['tile_offset_y']
    origin = np.array([x0, y0])
    df = df - origin
    return df.astype(int)


def coofy(spots, tile):
    spots_loc = localise_coords(spots, tile)
    x = spots_loc.x.values.astype(int)
    y = spots_loc.y.values.astype(int)

    coo_shape = np.array(tile['label_image'].shape)
    idx = spots.index.values + 1  # avoid having idx = 0
    coo = coo_matrix((idx, (y, x)), shape=coo_shape)

    return coo