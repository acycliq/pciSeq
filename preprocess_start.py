from src.preprocess.cell_merger import Stage
from src.preprocess.tile import Tile
import pandas as pd
from src.preprocess.spot_label import spot_label
from src.preprocess.cell_props import calc_props
import logging
import config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)
logger = logging.getLogger()

if __name__ == "__main__":
    cfg = config.PREPROCESS

    spots_full = pd.read_csv(cfg['spots'])
    tiles = Tile(cfg)
    stage = Stage(tiles, spots_full)

    stage.merge_cells()
    stage.post_merge([stage.tiles, stage.merge_register, stage.label_generator])
    stage.cell_props = calc_props(stage)

    stage.cell_props['cell_id'] = stage.assign_cell_id()

    logger.info('Spot labelling started')
    stage.spots = spot_label(stage)
    logger.info('Spot labelling done')

    # Save now the data on the filesystem
    stage.writer()

    print('Done!')
