"""
does the cell typing
"""
import os
from src.cell_call.main import VarBayes
import logging
import config

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

if __name__ == "__main__":
    case = 'MOUSE'  # 'MOUSE' or 'HUMAN'

    ini = getattr(config, case)
    out_dir = os.path.join(config.ROOT_DIR, 'data', 'cell_call_demo_data',  case, 'cell_type_output')

    # 1. run the cell calling algo
    varBayes = VarBayes(ini)
    cellData, geneData = varBayes.run()

    # 2. save the results
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cellData.to_json(os.path.join(out_dir, 'cellData_check.json'), orient='records')
    logger.info('Saved at %s' % (os.path.join(out_dir, 'cellData_check.json')))

    geneData.to_json(os.path.join(out_dir, 'geneData_check.json'), orient='records')
    logger.info('Saved at %s' % (os.path.join(out_dir, 'geneData_check.json')))


    print('Done')
