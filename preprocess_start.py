from src.preprocess.fov import Fov
from src.preprocess.cell_merger import Stage
from src.preprocess import utils
import os
import pandas as pd
# import data_manager.post as post
# import data_manager.result_splitter as rs
import logging
import config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

logger = logging.getLogger()



if __name__ == "__main__":
    # read the spots from the input file
    spots_all = pd.read_csv(os.path.join(config.ROOT_DIR, 'data', 'from_Matlab', 'SpotGlobal.csv'))

    fovs_across = config.PREPROCESSOR['FOVS_ACROSS']
    fovs_down = config.PREPROCESSOR['FOVS_DOWN']
    cfg = config.PREPROCESSOR

    filepath = os.path.join(config.ROOT_DIR, 'CellMap_left.mat')
    cellmap_chunks = utils.split_CellMap(filepath, config.PREPROCESSOR['fov_shape'][0], config.PREPROCESSOR['fov_shape'][1])
    fovs = Fov(cfg)
    stage = Stage(fovs, spots_all, cellmap_chunks)

    stage.merge_cells()
    stage.global_labels_par()
    stage.cell_props = stage.calc_props()
    stage.cell_props['cell_id'] = stage.assign_cell_id()

    logger.info('Spot labelling started')
    stage.spots = stage.assign_spot_parent()
    logger.info('Spot labelling done')

    # Save now the data on the filesystem
    stage.writer()

    # post.block_boundaries(fovs)
    # post.feature_collector(fovs, cfg)
    #
    # rs.split(stage.fovs, cfg)

    print('Done!')