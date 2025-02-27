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
    'InsideCellBonus': 0,

    # MisreadDensity: Expected number of misread spots. A dictionary contains user-defined values
    # for gene misread densities used in the analysis.
    # The process to determine the misread density for each gene is as follows:
    #
    # 1. Compute the misread density as the number of misreads per gene divided by the area of the image.
    # 2. Calculate the mean of these computed misread densities.
    # 3. If the computed mean is NaN (for example, due to insufficient data), the fallback value
    #    specified by the key 'default' in this configuration is used.
    # 4. Finally, update the computed densities with any gene-specific overrides provided here.
    #
    # The 'default' key thus sets a baseline misread density to ensure that every gene is assigned
    # a consistent value when the mean cannot be reliably determined.
    #
    # Example configuration:
    # {
    #     'default': 1e-06,  # Fallback misread density if the computed mean is NaN.
    #     'Plp1': 0.0001,    # User-defined misread density for gene 'Plp1'.
    #     # Additional gene-specific overrides can be added here.
    # }
    'MisreadDensity': 0.00001,

    # A pseudo-count representing the confidence in the cell centroid estimated
    # by the segmentation algorithm.
    #   - Think of `cell_centroid_prior` as if you had already seen this many imaginary
    #     data points, all centered at the segmentation-derived centroid, before
    #     incorporating real data.
    #   - A small value means you have little confidence in the segmentation result,
    #     so real data will strongly influence the estimated centroid.
    #   - A large value means you trust the segmentation strongly, so the estimated
    #     centroid will remain close to the initial segmentation result, even as new
    #     data is introduced.
    'cell_centroid_prior': 10,

    # A pseudo-count representing the confidence in the initial estimate of covariance,
    # which is modeled as a diagonal matrix where each diagonal element is the square
    # of the mean cell radius.
    #   - Imagine this prior as if you had already observed this many imaginary data points
    #     that reinforce your belief about the spread (covariance).
    #   - A small value means you are uncertain about the initial covariance estimate,
    #     allowing real data to significantly influence the updated covariance.
    #   - A large value means you strongly trust the initial covariance assumption, so
    #     real data will only cause gradual updates.
    'cell_cov_prior': 10,

    # Gene detection might come with irregularities due to technical errors. A small value is introduced
    # here to account for these errors. It is an additive factor, applied to the single cell expression
    # counts when the mean counts per class and per gene are calculated.
    # It is like a tiny safety cushion for gene counts and adds a tiny number to all our counts to 
    # help handle these small errors.
    # This is especially helpful when we see zero counts, as it prevents mathematical problems
    # when we're doing calculations with these numbers.
    'SpotReg': 0.1,

    # By default, only the 6 nearest cells will be considered as possible parent cells for any given spot.
    # There is also one extra 'super-neighbor', which is always a neighbor to the spots so we can assign
    # the misreads to. Could be seen as the background. Hence, by default the algorithm tries examines
    # whether any of the 3 nearest cells is a possible parent cell to a given cell or whether the spot is
    # a misread
    'nNeighbors': 6,

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

    # voxel_size: Physical size of voxels in each dimension [x, y, z].
    # Used to correct for anisotropic sampling in 3D data.
    # Example: For a microscope with:
    #   - xy resolution of 0.147 µm/pixel
    #   - z-step size of 0.9 µm
    # Use: [0.147, 0.147, 0.9]
    # Default: [1, 1, 1] (isotropic voxels)
    'voxel_size': [1, 1, 1],

    # exclude_planes: List of z-planes to exclude from pciSeq.
    # Useful for removing poor quality or out-of-focus planes.
    # Example: [0, 1, 2, 3, 79, 80, 81, 82, 83] excludes first four and last five planes
    # Note: After exclusion, z coordinates are adjusted relative to
    # the first remaining plane.
    # Default: None (use all planes)
    'exclude_planes': None,

    # remove_flat_cells: Controls removal of cells that appear in only one z-plane.
    # These single-plane cells are often artifacts from segmentation, especially in 3D data.
    # Note: Only relevant for 3D data (multiple z-planes).
    # For 2D data, this setting has no effect.
    'remove_flat_cells': True,

    # Runtime attribute (automatically set during execution)
    'is3D': None,

    # *******************************************************************************
    # Hyperparameters below come into action **ONLY** if single cell data are missing
    # *******************************************************************************
    'mean_gene_counts_per_class': 60,
    'mean_gene_counts_per_cell': 30,



}

