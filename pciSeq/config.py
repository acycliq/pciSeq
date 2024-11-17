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

    # MisreadDensity: Expected number of misread spots per unit area (2D) or volume (3D)
    # Can be set as either:
    #   - Scalar value: Same density applied everywhere
    #   - Dict: Different densities for specific genes
    #
    # Important: All coordinates are first normalized to pixel units using voxel_size.
    # Example: If voxel_size = [0.147, 0.147, 0.9]:
    # - Each Z step equals 6.12 pixels in normalized space (0.9/0.147)
    # - A volume of 100x100x10 raw pixels becomes 100x100x61.2 normalized pixels
    #
    # Option 1 - Scalar (same value for all genes):
    # 'MisreadDensity': 0.0001,    # 2D: 1 misread per 100x100 pixel area OR
    #                   0.000001   # 3D: 1 misread per 100x100x100 normalized pixel volume
    #                                      (~100x smaller values due to volume vs area)
    #
    # Option 2 - Dictionary (different values per gene):
    # For 2D data (per pixel area):
    # 'MisreadDensity': {
    #     'default': 1e-5,     # 2D: 1 misread per 100x100 pixel area. Used for any genes not explicitly listed
    #     'Vip': 1e-5,         # Clean signal, low background
    #     'Npy': 1e-4,         # Noisier gene with more background
    #     'Aldoc': 1e-6,       # Very clean signal
    # },
    #
    #
    # For 3D data (per pixel volume):
    # 'MisreadDensity': {
    #     'default': 1e-5,     # 3D: 100 misreads per 100x100x100 normalized pixel volume. Used for any genes not explicitly listed
    #     'Vip': 1e-5,         # Clean signal, low background
    #     'Npy': 1e-4,         # Noisier gene with more background
    #     'Aldoc': 1e-6,       # Very clean signal
    # },
    #
    # Note: Values might be similar between 2D and 3D if:
    # - Background noise is uniform throughout the tissue
    # - Number of spots scales proportionally with volume
    # - All coordinates are properly normalized using voxel_size
    #
    # Rules of thumb for setting values:
    # For 2D data (per pixel area):
    # - A value of 0.0001 means 1 misread per 100x100 pixel area (since 1/10000 = 0.0001)
    # - Default (0.00001) means 0.1 misread per 100x100 pixel area
    # - For high-quality data: use 0.000001 to 0.00001
    # - For noisier data: use 0.0001 to 0.001
    #
    # For 3D data (per pixel volume):
    # - Use ~100x smaller values than 2D due to volume vs area
    # - Default (0.00001) means 10 misread per 100x100x100 normalized pixel volume
    # - For high-quality data: use 0.0000001 to 0.000001
    # - For noisier data: use 0.00001 to 0.0001
    #
    # To calculate a custom value:
    # 2D Example:
    # - 10 background spots in 200x200 pixels
    # - Area = 200 * 200 = 40,000 pixels²
    # - Density = 10/40,000 = 0.00025 misreads per pixel²
    #
    # 3D Example with voxel_size=[0.147, 0.147, 0.9]:
    # - 10 background spots in 200x200x10 raw voxels
    # - Z dimension normalized: 10 * (0.9/0.147) = 61.2 pixels
    # - Normalized volume = 200 * 200 * 61.2 = 2,448,000 pixels³
    # - Density = 10/2,448,000 = 0.000004 misreads per pixel³
    #
    # Note: The example above assumes the same number of background
    # spots, hence the 3D density is naturally much lower because
    # it is divided by volume instead of area.
    # If background spots are uniformly distributed in the tissue,
    # the density values might end up similar between 2D and 3D!
    # Always adjust your density values based on whether you're
    # working with 2D or 3D data!
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
    'save_data': True,

    # Set here where the results will be saved. If default then they will be saved at your system's temp folder
    'output_path': 'default',

    # if true the viewer will be launched once convergence has been achieved
    'launch_viewer': False,

    'launch_diagnostics': False,

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

    # Runtime attribute (automatically set during execution)
    'is_redis_running': False,

    # *******************************************************************************
    # Hyperparameters below added for 3D
    # *******************************************************************************
    'voxel_size': [1, 1, 1],  # x, y, z

    'exclude_planes': None,

    # Runtime attribute (automatically set during execution)
    'is3D': None,

    # *******************************************************************************
    # Hyperparameters below come into action **ONLY** if single cell data are missing
    # *******************************************************************************
    'mean_gene_counts_per_class': 60,
    'mean_gene_counts_per_cell': 30,



}

