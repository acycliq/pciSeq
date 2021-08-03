"""
Standalone script to carve out the roi-specific cell by gene array
WARNING. THIS OPERATES ON THE NON ROTATED ROI
"""
import os
import glob
import numpy as np
import pandas as pd
import pciSeq.src.preprocess.vizgen.config as config
from pathlib import Path
from pciSeq.src.preprocess.vizgen.utils import transformation, is_inside
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def run(merfish_id, slice_id, region_id):
    logger.info('Doing %s, %s, %s' % (merfish_id, slice_id, region_id))
    cfg = config.get_config(merfish_id=merfish_id, slice_id=slice_id, region_id=region_id)

    # read the cell metadata
    cell_meta = pd.read_csv(cfg['cell_metadata'], index_col=0)

    # Transformation from micron to pixel coords
    tx, ty, _, _ = transformation(cfg)

    # get the pixel coordinates of the cell centroids
    cell_meta['center_x_px'] = tx(cell_meta.center_x.values)
    cell_meta['center_y_px'] = ty(cell_meta.center_y.values)

    # get now the roi
    roi = cfg['clip_poly']

    # find which cell centroid fall inside the roi
    points = cell_meta[['center_x_px', 'center_y_px']].values
    poly = np.array(roi)
    mask = is_inside(points, poly)

    # filter the full dataframe
    idx = cell_meta[mask].index.values
    df = pd.read_csv(cfg['cell_by_gene'], index_col=0).loc[idx]
    logger.info('Found %d cell inside the roi' % df.shape[0])

    dir_name = os.path.join(config.DROPBOX_URL, 'Dimitris folder', 'cell_by_gene', merfish_id, slice_id, region_id, 'roi_cell_by_gene')
    file_name = save_df_simple(df, dir_name)
    logger.info('cell by gene array saved to %s' % file_name)
    return df


def save_df_simple(df, dir_path):
    """
    Simple convenience function to save a dataframe as a tsv.
    It creates the containing folder if it doesnt exist
    :param df:
    :param dir_path:
    :return:
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    else:
        files = glob.glob(dir_path + '/*.*')
        for f in files:
            os.remove(f)
    filename = os.path.basename(dir_path).split('.')[0]
    file_name = os.path.join(dir_path, filename + '.csv')
    df.to_csv(file_name)
    return file_name

def get_slice_ids(merfish_id):
    return get_subdir(merfish_id)


def get_region_ids(merfish_id, slice_id):
    return get_subdir(merfish_id, slice_id)


def get_subdir(merfish_id, slice_id=''):
    fPath = os.path.join(config.DROPBOX_URL, merfish_id, slice_id)
    p = Path(fPath)
    out = [os.path.basename(f) for f in p.iterdir() if f.is_dir()]
    return out


if __name__ == "__main__":
    merfish_ids = ['MERFISH_F_E', 'MERFISH_F_F', 'MERFISH_M_C', 'MERFISH_M_Z']

    for merfish_id in merfish_ids:
        for slice_id in get_slice_ids(merfish_id):
            for region_id in get_region_ids(merfish_id, slice_id):
                # logger.info("\n Started slice %s, region %s" % (slice_id, region_id))
                run(merfish_id, slice_id, region_id)


    logger.info('Done!')

