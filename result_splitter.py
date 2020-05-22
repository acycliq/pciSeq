'''
postprocessing acript. Takes the cell type results, namely the json file that
describes which cell the spots belong to and splits it per fov
'''

import pandas as pd
from src.preprocess.fov import Fov
from src.preprocess.cell_merger import Stage,  get_dir
import os
import config
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
    )
logger = logging.getLogger()


def split(fovs, my_cfg):
    geneData = pd.read_json(my_cfg['CELL_TYPED_GENES'])
    cellData = pd.read_json(my_cfg['CELL_TYPED_CELLS'])

    # fovs_across = 20
    # fovs_down = 23
    #
    # fov_obj = Fov(fovs_across, fovs_down, cfg)

    for i, fov in enumerate(fovs):
        assert i == fov['fov_id']
        x_range = fov['fov_range']['x']
        y_range = fov['fov_range']['y']

        split_geneData(fov, geneData, x_range, y_range, my_cfg)
        split_cellData(fov, cellData, x_range, y_range, my_cfg)


def split_geneData(fov, geneData, x_range, y_range, my_cfg):
    mask = (geneData.x.values >= x_range[0]) & \
           (geneData.x.values < x_range[1]) & \
           (geneData.y.values >= y_range[0]) & \
           (geneData.y.values < y_range[1])
    df = geneData[mask]

    fov_id = fov['fov_id']
    fov_dir = get_dir(my_cfg, fov_id)
    full_path = os.path.join(fov_dir, 'cell_type_out', 'fov_' + str(fov_id) + '_Dapi_overlays.json')

    if not os.path.exists(os.path.dirname(full_path)):
        os.makedirs(os.path.dirname(full_path))

    # add an extra column and populate it with the fov_id
    df['fov_id'] = fov['fov_id']

    df.to_json(full_path, orient='records')
    logger.info('Fov_id: %d: Dapi overlays saved at %s' % (fov_id, full_path))


def split_cellData(fov, cellData, x_range, y_range, my_cfg):
    mask = (cellData.X.values >= x_range[0]) & \
           (cellData.X.values < x_range[1]) & \
           (cellData.Y.values >= y_range[0]) & \
           (cellData.Y.values < y_range[1])
    df = cellData[mask]

    # add an extra column and populate it with the fov_id
    df['fov_id'] = fov['fov_id']

    fov_id = fov['fov_id']
    fov_dir = get_dir(my_cfg, fov_id)
    full_path = os.path.join(fov_dir, 'cell_type_out', 'fov_' + str(fov_id) + '_iss.json')

    if not os.path.exists(os.path.dirname(full_path)):
        os.makedirs(os.path.dirname(full_path))

    df.to_json(full_path, orient='records')
    logger.info('Fov_id: %d: cell data saved at %s' % (fov_id, full_path))


if __name__ == "__main__":
    fovs_across = config.PREPROCESSOR['FOVS_ACROSS']
    fovs_down = config.PREPROCESSOR['FOVS_DOWN']
    cfg = config.PREPROCESSOR

    fovs = Fov(fovs_across, fovs_down, cfg)
    stage = Stage(fovs)

    split(stage.fovs, cfg)
