import os

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
print(ROOT_DIR)


PREPROCESS = {
    'spots': os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'iss', 'spots.csv'),
    'label_image': os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'segmentation', 'label_image.coo.npz'),

    # Optional setting. If this is set, then the label_image will be split into smaller arrays (tiles).
    # If it is set to [None, None] the tile dims will be overridden by the image dimensions
    'tile_size': [None, None],  # [width_px, height_px]

    # Target folder to save temp data from the preprocessing step
    'temp': os.path.join(ROOT_DIR, 'src', 'preprocess', 'temp')
}


MOUSE = {
    'CellCallTolerance': 0.02,
    'Inefficiency': 0.2,
    'InsideCellBonus': 2,
    'MisreadDensity': 0.00001,
    'SpotReg': 0.1,
    'nNeighbors': 3,
    'rGene': 20,
    'rSpot': 2,
    'max_iter': 100,
    'scRNAseq': os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'scRNA', 'scRNAseq.csv.gz'),
    # 'spotsFile': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse', 'spots_ext.csv'),  # Spot attributes, contains x,y coordinates for the spots and their gene names
    # 'expanded_cells':  os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse', 'expanded_cells', 'expanded_cells.csv'),
    'exclude_genes': []
}
