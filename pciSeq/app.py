"""
does the cell typing
"""
import os
from pciSeq.src.cell_call.main import VarBayes
from pciSeq.src.preprocess.spot_labels import stage_data
from pciSeq.src.viewer.utils import splitter_mb
import logging
# from pciSeq import config
from configparser import ConfigParser

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def cell_type(ini):
    # 1. run the cell calling algo
    varBayes = VarBayes(ini)
    cellData, geneData = varBayes.run()

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


def app(ini):
    cfg = ConfigParser()
    cfg.read(ini)

    # 1. prepare the data
    stage_data(cfg['PREPROCESS'])

    # 2. cell typing
    cell_type(cfg)  # 'MOUSE' or 'HUMAN'
    logger.info('Done')


if __name__ == "__main__":
    # Read config.ini file
    cfg = ConfigParser()
    cfg.read('./pciSeq/config.ini')

    # 1. prepare the data
    stage_data(cfg['PREPROCESS'])

    # 2. cell typing
    cell_type(cfg)  # 'MOUSE' or 'HUMAN'
    logger.info('Done')
