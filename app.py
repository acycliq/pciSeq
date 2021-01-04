"""
does the cell typing
"""
import os
from src.cell_call.main import VarBayes
from src.preprocess import stage_data
from src.viewer.utils import splitter_mb
import logging
import config

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def pciSeq(ini):
    # 1. run the cell calling algo
    varBayes = VarBayes(ini)
    cellData, geneData = varBayes.run()

    # 2. save the results
    out_dir = ini['out_dir']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cellData.to_csv(os.path.join(out_dir, 'cellData.tsv'), sep='\t', index=False)
    logger.info('Saved at %s' % (os.path.join(out_dir, 'cellData.tsv')))

    geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
    logger.info('Saved at %s' % (os.path.join(out_dir, 'geneData.tsv')))

    # Write to the disk as tsv of 99MB each
    splitter_mb(cellData, os.path.join(out_dir, 'cellData'), 99)
    splitter_mb(geneData, os.path.join(out_dir, 'geneData'), 99)


if __name__ == "__main__":
    # 1. prepare the data
    # stage_data.run(config.PREPROCESS)

    # 2. cell typing
    pciSeq(config.MOUSE)  # 'MOUSE' or 'HUMAN'
    logger.info('Done')
