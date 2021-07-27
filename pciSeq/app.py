import os
import pandas as pd
import numpy as np
from typing import Tuple
from pciSeq.src.cell_call.main import VarBayes
from pciSeq.src.preprocess.vizgen.stage_vizgen import stage_vizgen
from pciSeq.src.viewer.utils import splitter_mb
from pciSeq import config
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
DATASTORE_DIR = os.path.join('D:\\', 'rotated_dapi_map_tiles')

def fit(scRNAseq: pd.DataFrame, opts: dict = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point for pciSeq.

    Parameters
    ----------
    iss_spots : pandas.DataFrame
        Index:
            RangeIndex
        Columns:
            Name: Gene, dtype: string, The gene name
            Name: x, dtype: int64, X-axis coordinate of the spot
            Name: y, dtype: int64, Y-axis coordinate of the spot

    scRNAseq : pandas.DataFrame
        Index:
            The gene name
        Columns:
            The column headers are the cell classes and the data are uint32

    opts : dictionary (Optional)
        A dictionary to pass-in user-defined hyperparameter values. They override the default
        values as these are set by the config.py file. For example to exclude genes Npy and
        Vip you can define opts as:
            opts = {'exclude_genes': ['Npy', 'Vip']}
        and pass that dict to the fit function as the last argument

    Returns
    ------
    cellData : pandas.DataFrame
        Index:
            RangeIndex
        Columns:
            Name: Cell_Num, dtype: int64, The label of the cell
            Name: X, dtype: float64, X-axis coordinate of the cell centroid
            Name: Y, dtype: float64, Y-axis coordinate of the cell centroid
            Name: Genenames, dtype: Object, array-like of the genes assinged to the cell
            Name: CellGeneCount, dtype: Object,array-like of the corresponding gene counts
            Name: ClassName, dtype: Object, array-like of the genes probable classes for the cell
            Name: Prob, dtype: Object, array-like array-like of the corresponding cell class probabilities

    geneData : pandas.DataFrame
        Index:
            RangeIndex
        Columns:
            Name: Gene, dtype: string, The gene name
            Name: Gene_id, dtype: int64, The gene id, the position of the gene if all genes are sorted.
            Name: x, dtype: int64, X-axis coordinate of the spot
            Name: y, dtype: int64, Y-axis coordinate of the spot
            Name: neighbour, dtype: int64, the label of the cell which is more likely to 'raise' the spot. If zero then the spot is a misread.
            Name: neighbour_array, dtype: Object, array-like with the labels of the 4 nearest cell. The last is always the background and has label=0
            Name: neighbour_prob, dtype: Object, array-like with the prob the corresponding cell from neighbour_array has risen the spot.
    """

    # 1. get the hyperparameters
    cfg = init(opts)

    # 2. prepare the data
    logger.info(' Preprocessing data')
    _cells, cellBoundaries, _spots = stage_vizgen(cfg)
    logger.info(' Preprocessing finished')

    # 3. cell typing
    cellData, geneData = cell_type(_cells, _spots, scRNAseq, cfg)

    # 4. save to filesystem
    save_data = True
    if save_data:
        write_data(cellData, geneData, cellBoundaries, cfg)

    logger.info(' Done')
    return cellData, geneData


def cell_type(_cells, _spots, scRNAseq, ini):
    varBayes = VarBayes(_cells, _spots, scRNAseq, ini)

    logger.info(' Start cell typing')
    cellData, geneData = varBayes.run()
    return cellData, geneData


def write_data(cellData, geneData, cellBoundaries, ini):
    # out_dir = ini['out_dir']
    out_dir = os.path.join(DATASTORE_DIR, ini['slice_id'], ini['region_id'], 'celltyped')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cellData.to_csv(os.path.join(out_dir, 'cellData.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'cellData.tsv')))

    geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'geneData.tsv')))

    cellBoundaries.to_csv(os.path.join(out_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
    logger.info(' Saved at %s' % (os.path.join(out_dir, 'cellBoundaries.tsv')))

    # Write to the disk as tsv of 99MB each
    splitter_mb(cellData, os.path.join(out_dir, 'cellData'), 99)
    splitter_mb(geneData, os.path.join(out_dir, 'geneData'), 99)
    splitter_mb(cellBoundaries, os.path.join(out_dir, 'cellBoundaries'), 99)


def init(opts):
    """
    Reads the opts dict and if not None, it will override the default parameter value by
    the value that the dictionary key points to.
    If opts is None, then the default values as these specified in the config.py file
    are used without any change.
    """
    cfg = config.DEFAULT
    if opts is not None:
        default_items = set(cfg.keys())
        user_items = set(opts.keys())
        assert user_items.issubset(default_items), ('Options passed-in should be a dict with keys: %s ' % default_items)
        for item in opts.items():
            if isinstance(item[1], str):
                val = item[1]
            elif isinstance(item[1], (int, float, list)) or isinstance(item[1](1), np.floating):
                val = item[1]
            else:
                raise TypeError("Only integers, floats and lists and strings are allowed")
            cfg[item[0]] = val
            logger.info(' %s is set to %s' % (item[0], cfg[item[0]]))
    return cfg


if __name__ == "__main__":
    slice_id = "MsBrain_Eg3_VS6_JH_V6_05-01-2021"
    region_id = "region_1"

    # read some demo data
    _scRNAseq = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'vizgen', 'scRNA', 'scRNAseq.csv.gz'),
                            header=None, index_col=0, compression='gzip', dtype=object)
    _scRNAseq = _scRNAseq.rename(columns=_scRNAseq.iloc[0], copy=False).iloc[1:]
    _scRNAseq = _scRNAseq.astype(np.float).astype(np.uint32)

    # pass now the paths to the preprocessed vizgen data
    opts = {
        'slice_id': slice_id,
        'region_id': region_id,
        'cellBoundaries_file': os.path.join(DATASTORE_DIR, slice_id, region_id, 'cellBoundaries', 'cellBoundaries.tsv'),
        'cell_props_file': os.path.join(DATASTORE_DIR, slice_id, region_id, 'cell_props', 'cell_props.tsv'),
        'spots_file': os.path.join(DATASTORE_DIR, slice_id, region_id, 'labelled_spots', 'labelled_spots.tsv'),
        'clip_poly_file': os.path.join(DATASTORE_DIR, slice_id, region_id, 'roi', 'roi_rotated.csv'),
    }
    fit(_scRNAseq, opts)
