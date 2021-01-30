"""
does the cell typing
"""
import os
import pandas as pd
import numpy as np
from pciSeq.src.cell_call.main import VarBayes
from pciSeq.src.preprocess.spot_labels import stage_data
from pciSeq.src.viewer.utils import splitter_mb
import logging
from scipy.sparse import coo_matrix, save_npz, load_npz
from configparser import ConfigParser

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def cell_type(_cells, _spots, scRNAseq, ini):
    # 1. run the cell calling algo
    varBayes = VarBayes(_cells, _spots, scRNAseq, ini)
    cellData, geneData = varBayes.run()

    save_data = False
    if save_data:
        # 2. save the results
        out_dir = ini['PCISEQ']['out_dir']
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        cellData.to_csv(os.path.join(out_dir, 'cellData.tsv'), sep='\t', index=False)
        logger.info('Saved at %s' % (os.path.join(out_dir, 'cellData.tsv')))

        geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
        logger.info('Saved at %s' % (os.path.join(out_dir, 'geneData.tsv')))

        # Write to the disk as tsv of 99MB each
        splitter_mb(cellData, os.path.join(out_dir, 'cellData'), 99)
        splitter_mb(geneData, os.path.join(out_dir, 'geneData'), 99)
    return cellData, geneData


def app(iss_spots, coo, scRNAseq, cfg):
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

    coo : scipy.sparse.coo_matrix
        A label image array as a coo_matrix datatype. The label denote
        which cell the corresponding pixel 'belongs' to. If label is
        zero, the pixel is on the background

    scRNAseq : pandas.DataFrame
        Index:
            The gene name
        Columns:
            The column headers are the cell classes and the type of their values is np.uint32

    cfg : configParser
        A configParser object which has parsed the config.ini to read the hyperparameters.

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
            Name: neighbour, dtype: int64, the label of the cell which is more likely to 'raise' the spot
            Name: neighbour_array, dtype: Object, array-like with the labels of the 4 nearest cell. The last always the background and has label=0
            Name: neighbour_prob, dtype: Object, array-like with the prob the corresponding cell from neighbour_array has risen the spot.
    """

    # 1. prepare the data
    logger.info('Preprocessing data')
    _cells, _cell_boundaries, _spots = stage_data(iss_spots, coo)

    # 2. cell typing
    logger.info('Start cell typing')
    cellData, geneData = cell_type(_cells, _spots, scRNAseq, cfg)
    logger.info('Done')
    return cellData, geneData


if __name__ == "__main__":
    ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

    # read config.ini file
    cfg = ConfigParser()
    cfg.read(os.path.join(ROOT_DIR, 'config.ini'))

    # read some demo data
    _iss_spots = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'iss', 'spots.csv'))
    _coo = load_npz(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'segmentation', 'label_image.coo.npz'))

    _scRNAseq = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'scRNA', 'scRNAseq.csv.gz'),
                            header=None, index_col=0, compression='gzip', dtype=object)
    _scRNAseq = _scRNAseq.rename(columns=_scRNAseq.iloc[0], copy=False).iloc[1:]
    _scRNAseq = _scRNAseq.astype(np.float).astype(np.uint32)

    # main task
    app(_iss_spots, _coo, _scRNAseq, cfg)

