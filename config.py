import os
import platform

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
print(ROOT_DIR)


PREPROCESSOR = {
    'MATLAB_SPOTS': os.path.join(ROOT_DIR, 'data', 'from_Matlab', 'split', 'spots'),
    'FOV_ROOT': os.path.join(ROOT_DIR, 'data', 'fov'),
    'fov_size': 2000,   # implies that each fov is square with side length 2000px
    'FOVS_ACROSS': 11,
    'FOVS_DOWN': 14,
    # 'fovs_across': 11,  # This is not an input, it is calculated in '''split_image'''. I put it here for convenience
    # 'fovs_down': 14,    # This is not an input, it is calculated in '''split_image'''. I put it here for convenience
}


MOUSE = {
    'CellCallMaxIter': 100,
    'CellCallTolerance': 0.02,
    'Inefficiency': 0.2,
    'InsideCellBonus': 2,
    'MisreadDensity': 0.00001,
    'SpotReg': 0.1,
    'nNeighbors': 3,
    'rGene': 20,
    'rSpot': 2,
    'max_iter': 100,
    'scRNAseq': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse', 'scRNA', 'scRNAseq.csv.gz'),
    'saFile': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse', 'spots_ext.csv'),  # Spot attributes, contains x,y coordinates for the spots and their gene names
    'expanded_cells':  os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse', 'expanded_cells', 'expanded_cells.csv'),
    'exclude_genes': [],
    'drop_nan': False,  # That is only a temp solution. I will be removing that when everything is set correctly
}

MOUSE_FULL_CORONAL = {
    'CellCallMaxIter': 100,
    'CellCallTolerance': 0.02,
    'Inefficiency': 0.2,
    'InsideCellBonus': 2,
    'MisreadDensity': 0.00001,
    'SpotReg': 0.1,
    'nNeighbors': 3,
    'rGene': 20,
    'rSpot': 2,
    'max_iter': 100,
    'scRNAseq': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse_full_coronal', 'cell_type_input', 'scRNAseq.csv.gz'),
    'saFile': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse_full_coronal', 'cell_type_input', 'spots.csv'),  # Spot attributes, contains x,y coordinates for the spots and their gene names
    'expanded_cells':  os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse_full_coronal', 'cell_type_input', 'expanded_cells.csv'),
    'exclude_genes': [],
    'drop_nan': False,  # That is only a temp solution. I will be removing that when everything is set correctly
}