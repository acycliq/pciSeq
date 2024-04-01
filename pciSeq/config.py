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
    'rGene': 20,
    'Inefficiency': 0.2,

    # If a spot is inside the cell boundaries this bonus will give the likelihood an extra boost
    # in order to make the spot more probable to get assigned to the cell than another spot positioned
    # outside the cell boundaries
    'InsideCellBonus': 2,

    # To account for spots far from the some a uniform distribution is introduced to describe those misreads.
    # By default this uniform distribution has a density of 1e-5 misreads per pixel.
    'MisreadDensity': 0.00001,

    # Gene detection might come with irregularities due to technical errors. A small value is introduced
    # here to account for these errors. It is an additive factor, applied to the single cell expression
    # counts when the mean counts per class and per gene are calculated.
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
    # Hyperparameters below come into action **ONLY** if single cell data are missing
    # *******************************************************************************
    'mean_gene_counts_per_class': 60,
    'mean_gene_counts_per_cell': 15,

}

