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

