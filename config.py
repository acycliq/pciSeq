import os
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)
logger = logging.getLogger()

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
logger.info('Root dir is: %s' % ROOT_DIR)

# Settings to prepare the data. The label image and spots are parsed and if a spot
# lies within the cell boundaries then the corresponding cell id is recorded.
PREPROCESS = {
    'spots': os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'iss', 'spots.csv'),
    'label_image': os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'segmentation', 'label_image.coo.npz'),

    # Optional setting. If this is set, then the label_image will be split into smaller arrays (tiles).
    # If it is set to [None, None] the tile dims will be set equal to the image dimensions
    'tile_size': [None, None],  # [width_px, height_px]

    # Target folder to save temp data from the preprocessing step
    'temp': os.path.join(ROOT_DIR, 'out', 'temp')
}


MOUSE = {
    'exclude_genes': [],  # list of genes to be excluded from the cell type algo
    'out_dir': os.path.join(ROOT_DIR, 'out'),  # folder to save the results
    'scRNAseq': os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'scRNA', 'scRNAseq.csv.gz'),  # Single cell data

    # hyperparameters for the pciSeq method
    'CellCallTolerance': 0.02,
    'Inefficiency': 0.2,
    'InsideCellBonus': 2,
    'MisreadDensity': 0.00001,
    'SpotReg': 0.1,
    'nNeighbors': 3,
    'rGene': 20,
    'rSpot': 2,
    'max_iter': 100
}
