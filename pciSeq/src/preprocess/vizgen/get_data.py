"""
Convenience functions only to collect the data needed for the draft viewer
It reads the vizgen raw data and creates geneData, cellData and cellBoundaries needed
to visualise the data. It DOES NOT DO ANY CELLTYPING

Note:
I can some cases where two distinct cells have the same cell key. I do not know why, I havent looked into it
Look for example the cell keys:
190964149170372246785994222967545698003
83460, 83461
or
192975879491480759827730747397591978672
83462, 83463
"""

import os
import glob
import numpy as np
import pandas as pd
import pciSeq.src.preprocess.vizgen.config as config
import logging
from pciSeq.src.preprocess.vizgen.utils import transformation, splitter_mb, save_df, save_df_simple, rotate_data, is_inside_sm_parallel
from pciSeq.src.preprocess.vizgen.cellBorders import cell_boundaries_px_par
from pciSeq.src.preprocess.vizgen import label_spots

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)



def clip_data(df, cfg):
    """
    Keeps spots inside the area defined by predetermined polygon only
    :param df:
    :param cfg:
    :return:
    """
    if cfg['clip_poly']:
        logger.info('Found clipping poly. Keeping data inside %s' % cfg['clip_poly'])
        coords = cfg['clip_poly']
        if not coords[0] == coords[-1]:
            logger.info('Closing the polygon')
            coords.append(coords[0])

        points = df[['x', 'y']].values
        poly = np.array(coords)
        mask = is_inside_sm_parallel(points, poly)
        logger.info('Collected %d points inside ROI' % mask.sum())
        return df[mask]
    else:
        return df


def get_gene_data(cfg):
    # Transformation from micron to pixel coords
    tx, ty, _, _ = transformation(cfg)

    # read the spots from the raw files
    spots_path = cfg['detected_transcripts']
    chunks = pd.read_csv(spots_path, chunksize=100000)
    data = pd.concat(chunks)

    # data_z3 = data[data.global_z == 3]
    data_z3 = data
    logger.info('Getting the unique gene ids')
    [unique_gene_names, gene_id, counts_per_gene] = np.unique(data_z3.gene.values, return_inverse=True,
                                                              return_counts=True)
    logger.info('ok. done')

    data_z3 = data_z3.assign(global_x_px=tx(data_z3.global_x.values).astype(np.int32))
    data_z3 = data_z3.assign(global_y_px=ty(data_z3.global_y.values).astype(np.int32))
    data_z3['Gene_id'] = gene_id
    data_z3 = data_z3.sort_values(by=['global_x_px', 'global_y_px'])

    data_z3 = data_z3[['gene', 'global_x_px', 'global_y_px', 'global_z', 'Gene_id', 'fov']].rename(
        columns={'gene': "Gene", 'global_x_px': "x", 'global_y_px': "y", 'global_z': 'z'})

    data_z3['neighbour'] = np.ones(len(gene_id)).astype(np.int32)
    data_z3['neighbour_array'] = [[1] for i in range(len(gene_id))]
    data_z3['neighbour_prob'] = [[1.0] for i in range(len(gene_id))]

    # data_z3 = data_z3.drop([['Unnamed: 0', 'Unnamed: 0.1']], axis=1)

    logger.info('Data size: %d, %d' % (data_z3.shape[0], data_z3.shape[1]))

    logger.info('Sorting by x, y')
    data_z3 = data_z3.sort_values(by=['x', 'y'])
    logger.info('ok. done')
    # splitter_mb(data_z3, 'geneData', 99)
    data_z3 = data_z3.astype({'x': np.int32,
                              'y': np.int32,
                              'z': np.int32,
                              'fov': np.uint32,
                              'neighbour': np.int32})
    return data_z3


def roi_spots(cfg):
    """
    Collects the spots inside the ROI and then rotates them
    :param cfg:
    :param out_path:
    :return:
    """

    # 1. First fetch the data
    geneData = get_gene_data(cfg)

    # 2. Keep only spots inside the ROI
    logger.info('Keeping only the ROI spots')
    geneData = clip_data(geneData, cfg)

    # 3. Rotate the spots
    logger.info('Rotating the transcript data by %d degrees' % cfg['rotation'][0])
    rot = rotate_data(geneData[['x', 'y']].values.copy(), cfg)
    geneData.x = rot[:, 0].astype(np.int32)
    geneData.y = rot[:, 1].astype(np.int32)

    # 4. Save the spots to the disk
    # splitter_mb(geneData, os.path.join(out_path, 'geneData'), 99)
    save_df(geneData, os.path.join(cfg['target_dir'], 'geneData'))
    logger.info('Gene data saved at: %s' % os.path.join(cfg['target_dir'], 'geneData'))

    # 5. Finally save the rotated roi
    rotated_roi = rotate_roi(cfg)
    rotated_roi = pd.DataFrame(rotated_roi, columns=['x', 'y'])
    str_path = os.path.join(cfg['target_dir'], 'roi')
    if not os.path.exists(str_path):
        os.makedirs(str_path)
    else:
        files = glob.glob(str_path + '/*.*')
        for f in files:
            os.remove(f)
    rotated_roi.to_csv(os.path.join(str_path, 'roi_rotated.csv'), index=False)


def cell_boundaries(cfg):
    """
    Collects the rotated outer ring of the cells boundaries. It get all cells not only those
    within the ROI
    :param cfg:
    :param out_path:
    :return:
    """

    # collect all cell boundaries
    px_boundaries = cell_boundaries_px_par(cfg)

    # make a dataframe with the cell centroids and cell area
    cell_props = px_boundaries[['cell_key',	'cell_label', 'x', 'y',	'cell_area']]
    save_df_simple(cell_props, os.path.join(cfg['target_dir'], 'cell_props'))

    # # make a dataframe with only the boundaries
    boundaries = px_boundaries[['cell_label', 'cell_boundaries']]
    save_df_simple(boundaries, os.path.join(cfg['target_dir'], 'cellBoundaries'))
    logger.info('cell data saved at: %s' % os.path.join(cfg['target_dir']))


def rotate_roi(cfg):
    """
    rotates the clipping polygon
    :param cfg:
    :return:
    """
    return rotate_data(cfg['clip_poly'], cfg).astype(np.int32)


def run(merfish_id, slice_id, region_id):
    cfg = config.get_config(merfish_id=merfish_id, slice_id=slice_id, region_id=region_id)

    # 1. get the ROI spots (rotated)
    roi_spots(cfg)

    # 2. get all the cell boundaries (rotated)
    cell_boundaries(cfg)

    # 3. label now the spots
    label_spots.run(cfg)


