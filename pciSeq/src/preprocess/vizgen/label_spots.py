"""
Code to label the spots. It checks if a spot lies within the cell boundaries and if so then it takes note by
adding the cell key (the vizgen cell id) to the 'inside_cell_key' column of the spots dataframe. A cell key
of None means that the spot is on the background (it is not inside any cell)

Note the the code that uses parallelism may give slightly different spots labels compared to the non-parallel code.
That happens on spots that are exactly on the boundary between two cells.  Depending on how the processes that
execute the task happen in time and which cell was executed last the cell label may got overwritten from one adjacent
cell to another, and may be end up to a different (but the neighbouring) cell
"""
import pandas as pd
import numpy as np
import glob
from os import listdir
from os.path import isfile, join
from natsort import natsorted
import pciSeq.src.preprocess.vizgen.config as config
from pciSeq.src.preprocess.vizgen.utils import is_inside_sm_parallel
from multiprocessing import Pool, cpu_count
from functools import partial
import os
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def get_spots(cfg):
    folder = os.path.join(cfg['target_dir'], 'geneData')
    tsv_files = [f for f in listdir(folder) if isfile(join(folder, f)) and os.path.splitext(f)[1] == '.tsv']
    tsv_files = natsorted(tsv_files)
    df_list = []
    for tsv_file in tsv_files:
        fName = os.path.join(folder, tsv_file)
        df = pd.read_csv(fName, sep='\t')
        df_list.append(df)
    out = pd.concat(df_list, ignore_index=True)
    return out


def worker(cell_meta, cellBoundaries, spots, all_keys, fov):
    logger.info('Started fov: %d' % fov)
    # get the cell keys of the cells in the fov
    fov_cells_keys = cell_meta[cell_meta.fov == fov]

    # some of the cell keys may not be in the cell boundaries dataframe (probably because I looked on in z3 to capture the boundaries)
    missing = set(fov_cells_keys.index.values) - all_keys

    # remove the missing
    idx = set(fov_cells_keys.index.values) - missing

    # get the cells that are captured by this fov
    cells_in_fov = cellBoundaries.loc[list(idx)]

    spots_in_fov = spots[spots.fov == fov]
    points = spots_in_fov[['x', 'y']].values
    for cell_key, row in cells_in_fov.iterrows():
        outer_ring = row['cell_boundaries']
        if not outer_ring[-1] == outer_ring[0]:
            outer_ring.append(outer_ring[0])
        # find if a spots lies within the poly
        poly = np.array(outer_ring)
        mask = is_inside_sm_parallel(points, poly)

        # get the index of those spots that fall inside the poly
        spots_idx = spots_in_fov.index[mask]

        # backfill the inside_cell_key column with the corresponding key
        spots.loc[spots_idx, ['inside_cell_key']] = cell_key
    return spots.inside_cell_key.values


def run(cfg):
    # 1. read the rotated spots
    spots = get_spots(cfg)
    spots['inside_cell_key'] = np.empty(spots.shape[0], dtype=object)
    logger.info(spots)

    # get the cell metadata
    cell_meta = pd.read_csv(cfg['cell_metadata'], index_col=0)

    # get the (rotated) boundaries
    cellBoundaries = pd.read_csv(os.path.join(cfg['target_dir'], 'cellBoundaries', 'cellBoundaries.tsv'), sep='\t')
    cellProps = pd.read_csv(os.path.join(cfg['target_dir'], 'cell_props', 'cell_props.tsv'), sep='\t')

    # the boundaries appear to be strings. Convert them to a list of tuples
    _cell_boundaries = [eval(d) for d in cellBoundaries.cell_boundaries]
    for i, x in enumerate(_cell_boundaries):
        _cell_boundaries[i] = [tuple(d) for d in x]
    cellBoundaries.cell_boundaries = _cell_boundaries

    # add a cell key column to cellBoundaries dataframe
    assert np.all(cellProps.cell_label == cellBoundaries.cell_label), 'cellBoundaries and cell_props are not aligned'
    cellBoundaries['cell_key'] = cellProps.cell_key
    cellBoundaries = cellBoundaries.set_index('cell_key')
    all_keys = set(cellBoundaries.index.values)

    # 2. loop over the fovs
    fovs = np.unique(spots.fov.values)
    logger.info('start spot labelling')

    processes = cpu_count()
    print(f'Using {processes} processes')
    pool = Pool(processes=processes)
    res = pool.map(partial(worker, cell_meta, cellBoundaries, spots, all_keys),  fovs)
    pool.close()
    pool.join()

    # the results from the multiprocessing step is a list of list. Complile all resultes into one list
    out = np.empty(spots.shape[0], dtype=object)
    for lst in res:
        idx = [i for i, val in enumerate(lst) if val]
        out[idx] = lst[idx]

    # populate the inside_cell_key column
    spots['inside_cell_key'] = out

    # save now to the disk
    path_str = os.path.join(cfg['target_dir'], 'labelled_spots')
    if not os.path.exists(path_str):
        os.makedirs(path_str)
    else:
        files = glob.glob(path_str + '/*.*')
        for f in files:
            os.remove(f)
    spots[['Gene', 'x', 'y', 'z', 'Gene_id', 'inside_cell_key']].to_csv(os.path.join(path_str, 'labelled_spots.tsv'), sep='\t', index=False)
    logger.info('Labelled spots save at %s' % path_str)



if __name__ == "__main__":
    merfish_id = 'MERFISH_F_E'
    slice_id = 'MsBrain_Eg3_VS6_JH_V6_05-01-2021'
    region_id = 'region_1'
    cfg = config.get_config(merfish_id=merfish_id, slice_id=slice_id, region_id=region_id)
    run(cfg)
    logger.info('Done!')



    # for fov in fovs[:10]:
    #     logger.info('Started fov: %d' % fov)
    #     # get the cell keys of the cells in the fov
    #     fov_cells_keys = cell_meta[cell_meta.fov == fov]
    #
    #     # some of the cell keys may not be in the cell boundaries dataframe (probably because I looked on in z3 to capture the boundaries)
    #     missing = set(fov_cells_keys.index.values) - all_keys
    #
    #     # remove the missing
    #     idx = set(fov_cells_keys.index.values) - missing
    #
    #     # get the cells that are captured by this fov
    #     cells_in_fov = cellBoundaries.loc[list(idx)]
    #
    #     spots_in_fov = spots[spots.fov == fov]
    #     points = spots_in_fov[['x', 'y']].values
    #     for cell_key, row in cells_in_fov.iterrows():
    #         outer_ring = row['cell_boundaries']
    #         if not outer_ring[-1] == outer_ring[0]:
    #             outer_ring.append(outer_ring[0])
    #         # find if a spots lies within the poly
    #         poly = np.array(outer_ring)
    #         mask = is_inside_sm_parallel(points, poly)
    #
    #         # get the index of those spots that fall inside the poly
    #         spots_idx = spots_in_fov.index[mask]
    #
    #         # backfill the inside_cell_key column with the corresponding key
    #         spots.loc[spots_idx, ['inside_cell_key']] = cell_key
    #
    # spots.to_csv('workbench_par.tsv', sep='\t', index=False)
    # logger.info('finished spot labelling')

