"""
hyperparameters for the pciSeq method
"""
import numpy as np

DEFAULT = {

    # list of genes to be excluded during cell-typing, e.g ['Aldoc', 'Id2'] to exclude all spots from Aldoc and Id2
    'exclude_genes': [],

    # Maximum number of loops allowed for the Variational Bayes to run
    'max_iter': 1000,

    # Convergence achieved if assignment probabilities between two successive loops is less than the tolerance
    'CellCallTolerance': 0.02,

    # A gamma distribution expresses the efficiency of the in-situ sequencing for each gene. It tries to capture
    # the ratio of the observed over the theoretical counts for a given gene. rGene controls the variance and
    # Inefficiency is the average of this assumed Gamma distribution
    #
    # Example: If you expect 100 RNA molecules of a gene:
    # - With Inefficiency=0.2, you'll detect about 20 on average (20% detection rate)
    # - rGene=20 means this efficiency is fairly consistent between genes
    # - A lower rGene would mean more variance (e.g., some genes at 5% efficiency, others at 35%)
    #
    # This helps account for systematic differences in detection efficiency between genes
    # when making cell type assignments.
    'rGene': 20,
    'Inefficiency': 0.2,

    # If a spot is inside the cell boundaries this bonus will give the likelihood an extra boost
    # in order to make the spot more probable to get assigned to the cell than another spot positioned
    # outside the cell boundaries
    'InsideCellBonus': 2,

    # MisreadDensity accounts for RNA spots that appear far from any cell (background noise).
    # It represents the expected number of misread spots per pixel in your image.
    #
    # Rules of thumb for setting MisreadDensity:
    # - Default (0.00001) means 1 misread per 100x100 pixel area
    # - For high-quality data with clean background: use 0.000001 to 0.00001
    # - For noisier data with more background spots: use 0.0001 to 0.001
    # - If your data has many spots far from cells: increase this value
    # - If most spots are clearly associated with cells: decrease this value
    #
    # To calculate a custom value:
    # 1. Count the number of spots you think are background noise in a region
    # 2. Divide by the area (in pixels) of that region
    # Example: 10 background spots in a 200x200 pixel region = 10/(200*200) = 0.00025
    'MisreadDensity': 0.00001,

    # Gene detection might come with irregularities due to technical errors. A small value is introduced
    # here to account for these errors. It is an additive factor, applied to the single cell expression
    # counts when the mean counts per class and per gene are calculated.
    # It is like a tiny safety cushion for gene counts and adds a tiny number to all our counts to 
    # help handle these small errors.
    # This is especially helpful when we see zero counts, as it prevents mathematical problems
    # when we're doing calculations with these numbers.
    'SpotReg': 0.1,

    # By default only the 3 nearest cells will be considered as possible parent cells for any given spot.
    # There is also one extra 'super-neighbor', which is always a neighbor to the spots so we can assign
    # the misreads to. Could be seen as the background. Hence, by default the algorithm tries examines
    # whether any of the 3 nearest cells is a possible parent cell to a given cell or whether the spot is
    # a misread
    'nNeighbors': 3,

    # A gamma distributed variate from Gamma(rSpot, 1) is applied to the mean expression, hence the counts
    # are distributed according to a Negative Binomial distribution.
    # The value for rSpot will control the variance/dispersion of the counts
    # rSpot controls how much variation we expect to see in gene counts between cells of the same type.
    # It's used in a Negative Binomial distribution which models gene expression.
    #
    # Rules of thumb for setting rSpot:
    # - Default (2) is good for typical single-cell RNA data
    # - Lower values (0.5-1) mean high variability between cells
    #   → Use when you expect cells of the same type to show very different expression levels
    #   → Good for genes that tend to burst in expression
    # - Higher values (3-5) mean less variability between cells
    #   → Use when you expect cells of the same type to have similar expression levels
    #   → Good for housekeeping genes or very stable markers
    #
    # Examples:
    # rSpot = 0.5: Counts might vary a lot (e.g., [0,5,20,100] for same cell type)
    # rSpot = 2.0: Moderate variation (e.g., [10,15,20,25] for same cell type)
    # rSpot = 5.0: More consistent counts (e.g., [17,18,19,21] for same cell type)
    'rSpot': 2,

    # Boolean, if True the output will be saved as tsv files in a folder named 'pciSeq' in your system's temp dir.
    'save_data': False,

    # Set here where the results will be saved. If default then they will be saved at your system's temp folder
    'output_path': ['default'],

    # if true the viewer will be launched once convergence has been achieved
    'launch_viewer': False,

    'launch_diagnostics': True,

    # Initialise this to False, the correct value is set internally by the code itself
    'is_redis_running': False,

    # cell radius. If None then pciSeq will calc that as the mean radius across all cells.
    # Otherwise it will use the value provided below
    'cell_radius': None,

    # cell type prior: The prior distribution on the classes. It expresses the view on
    # how likely each class is to occur a-priori. It can be either 'uniform' or 'weighted'
    # 'uniform' means that the Zero class gets 50% and the remaining 50% is equally split
    # on the cell classes.
    # 'weighted' means that the cell type which is more likely to occur will be given more
    # weight. These weights are calculated dynamically within the algorithm based on
    # a Dirichlet distribution assumption.
    'cell_type_prior': 'uniform',

    # *******************************************************************************
    # Hyperparameters below added for 3D
    # *******************************************************************************
    'voxel_size': [1, 1, 1],  # x, y, z

    'exclude_planes': None,

    # this will be set automatically by the code
    'is3D': None,

    # *******************************************************************************
    # Hyperparameters below come into action **ONLY** if single cell data are missing
    # *******************************************************************************
    'mean_gene_counts_per_class': 60,
    'mean_gene_counts_per_cell': 30,



}

