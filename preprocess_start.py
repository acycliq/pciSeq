from src.preprocess.fov import Fov
from src.preprocess.cell_merger import Stage
from src.preprocess import utils
from src.preprocess.tile import Tile
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
    cfg = config.PREPROCESSOR

    spots_full = pd.read_csv(cfg['spots_full'])
    tiles = Tile(cfg)
    stage = Stage(tiles, spots_full)

    stage.merge_cells()
    stage.post_merge([stage.tiles, stage.merge_register, stage.label_generator])
    stage.cell_props = stage.calc_props()

    stage.cell_props['cell_id'] = stage.assign_cell_id()

    # Save now the data on the filesystem
    stage.writer()

    # post.block_boundaries(fovs)
    # post.feature_collector(fovs, cfg)
    #
    # rs.split(stage.fovs, cfg)

    print('Done!')