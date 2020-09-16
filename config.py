import os
import platform

ROOT_DIR = os.path.dirname(os.path.realpath(__file__))
print(ROOT_DIR)


PREPROCESSOR = {
    'fov_shape': [2000, 2000],  # implies that each fov has length 2000px and height 2000px
    'fovs_across': 11,
    'fovs_down': 14,
    'spots_full': os.path.join(ROOT_DIR, 'data', 'from_Matlab', 'SpotGlobal.csv'),
    'cellmap_full': os.path.join(ROOT_DIR, 'CellMap_left.mat'),

     # Output. Save here the results from the algorithm (I should move that somewhere else...)
    'CELL_TYPED_GENES': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse_full_coronal', 'cell_type_output', 'geneData.json'),
    'CELL_TYPED_CELLS': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse_full_coronal', 'cell_type_output', 'cellData.json'),
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
    # Hyperparameters
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

    # Inputs to the cell calling algorithm
    'scRNAseq': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse_full_coronal', 'cell_type_input', 'scRNAseq.csv.gz'),
    'saFile': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse_full_coronal', 'cell_type_input', 'spots.csv'),  # Spot attributes, contains x,y coordinates for the spots and their gene names
    'expanded_cells':  os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'mouse_full_coronal', 'cell_type_input', 'expanded_cells.csv'),
    'exclude_genes': [], # Write here the genes you want to exclude from the cell calling algorithm
    'drop_nan': False,  # That is only a temp solution. I will be removing that when everything is set correctly, Well remove it now!!
}

DAVID = {
    # Hyperparameters
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

    # Inputs to the cell calling algorithm
    'scRNAseq': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'David', 'scRNAseq.csv.gz'),
    'saFile': os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'David', 'spots.csv'),  # Spot attributes, contains x,y coordinates for the spots and their gene names
    'expanded_cells':  os.path.join(ROOT_DIR, 'data', 'cell_call_demo_data', 'David', 'cells.csv'),
    'exclude_genes': [], # Write here the genes you want to exclude from the cell calling algorithm
    'drop_nan': False,  # That is only a temp solution. I will be removing that when everything is set correctly, Well remove it now!!
}
