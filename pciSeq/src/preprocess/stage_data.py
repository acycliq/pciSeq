import pandas as pd
from pciSeq.src.preprocess.cell_merger import Stage
from pciSeq.src.preprocess.tile import Tile
from pciSeq.src.preprocess.spot_label import spot_label
from pciSeq.src.preprocess.cell_props import calc_props
import logging

logger = logging.getLogger(__name__)


def run(cfg):
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
    cells_df, cellCoords_df, spots_df = stage.writer(cfg['temp'])

    return cells_df, cellCoords_df, spots_df


if __name__ == "__main__":
    cfg = config.PREPROCESS
    run(cfg)
    print('Done!')

