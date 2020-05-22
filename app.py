import json
import os
import pandas as pd
from flask import Flask, render_template
import webbrowser
import platform
import random
from threading import Timer
from src.cell_call.run import varBayes
import pyvips
import shutil
import logging
import config

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

if __name__ == "__main__":
    case = 'MOUSE_FULL_CORONAL'  # 'MOUSE' or 'HUMAN'
    use_cache = False  # <--- I think i will remove that...

    my_config = getattr(config, case)
    out_dir = os.path.join(config.ROOT_DIR, 'data', 'cell_call_demo_data',  case, 'cell_type_output')

    # 1. run the cell calling algo
    if use_cache:  # I said it above, I am saying it here too: I think i will remove that...
        try:
            cellData = pd.read_csv(os.path.join(out_dir, 'cellData.csv'))
            geneData = pd.read_csv(os.path.join(out_dir, 'geneData.csv'))
            logger.info('Cached data loaded from %s' % out_dir)
        except IOError:
            logger.info('Could not read cache')
            cellData, geneData = varBayes(my_config)
    else:
        cellData, geneData = varBayes(my_config)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cellData.to_json(os.path.join(out_dir, 'cellData.json'), orient='records')
    logger.info('Saved at %s' % (os.path.join(out_dir, 'cellData.json')))

    geneData.to_json(os.path.join(out_dir, 'geneData.json'), orient='records')
    logger.info('Saved at %s' % (os.path.join(out_dir, 'geneData.json')))


    print('Done')
