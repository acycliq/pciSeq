"""
does the cell typing
"""
import os
from src.cell_call.main import VarBayes
from preprocess_start import run
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

    cellData.to_json(os.path.join(out_dir, 'cellData.json'), orient='records')
    logger.info('Saved at %s' % (os.path.join(out_dir, 'cellData.json')))

    geneData.to_json(os.path.join(out_dir, 'geneData.json'), orient='records')
    logger.info('Saved at %s' % (os.path.join(out_dir, 'geneData.json')))


if __name__ == "__main__":
    case = 'MOUSE'  # 'MOUSE' or 'HUMAN'
    ini = getattr(config, case)

    run(config.PREPROCESS)
    pciSeq(ini)
    print('Done')
