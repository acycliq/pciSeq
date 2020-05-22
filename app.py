'''
does the cell typing
'''
import os
from src.cell_call.run import varBayes
import logging
import config

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

if __name__ == "__main__":
    case = 'MOUSE_FULL_CORONAL'  # 'MOUSE' or 'HUMAN'

    my_config = getattr(config, case)
    out_dir = os.path.join(config.ROOT_DIR, 'data', 'cell_call_demo_data',  case, 'cell_type_output')

    # 1. run the cell calling algo
    cellData, geneData = varBayes(my_config)

    # 2. save the results
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cellData.to_json(os.path.join(out_dir, 'cellData.json'), orient='records')
    logger.info('Saved at %s' % (os.path.join(out_dir, 'cellData.json')))

    geneData.to_json(os.path.join(out_dir, 'geneData.json'), orient='records')
    logger.info('Saved at %s' % (os.path.join(out_dir, 'geneData.json')))


    print('Done')
