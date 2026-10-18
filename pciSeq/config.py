"""
Settings for pciSeq.

Pass any of these to pciSeq.fit() in an opts dictionary. Anything you leave out
keeps the default below.
"""


DEFAULT = {


    # Genes to exclude from cell calling, for example ['Aldoc', 'Id2'].
    "exclude_genes": [],


    # Maximum number of iterations. The loop ends here if the CellCallTolerance
    # criterion has not been met.
    "max_iter": 1000,


    # Stopping threshold. Iteration stops when the maximum absolute change in the
    # spot-to-cell assignment probabilities between two consecutive iterations,
    # taken over all spots and candidate cells, is below this value.
    "CellCallTolerance": 0.02,


    # Global scale factor on the cell type definitions, and the first
    # hyperparameter to tune on a new dataset. Every mean expression
    # value is multiplied by it before cells are scored, so it sets the overall
    # number of spots the model expects per cell. At 0.2, a gene with mean
    # expression 100 in a cell type is expected to give about 20 spots in a cell
    # of that type.
    #
    # Higher values raise the expected counts, so cells with few spots are more
    # often assigned to the Zero class. Lower values have the opposite effect.
    #
    # Inefficiency is the first of four multiplicative factors on the expected
    # counts, from the broadest to the most specific:
    #
    # | factor       | scope                        | prior strength |
    # | ------------ | ---------------------------- | -------------- |
    # | Inefficiency | all genes in all cells       | fixed          |
    # | eta          | a given gene, all cells      | rGene          |
    # | theta        | a given cell, all genes      | rTheta         |
    # | gamma        | a given gene in a given cell | rSpot          |
    #
    # Inefficiency is set directly. eta, theta and gamma are estimated from the
    # data, each with a prior of mean 1 whose strength is set by the
    # corresponding hyperparameter.
    "Inefficiency": 0.2,


    # Prior strength for eta, the per-gene detection factor. eta is the ratio of
    # the gene's observed count to its expected count over all cells, with a
    # Gamma(rGene, rGene) prior of mean 1. The estimate is
    # `eta = (observed + rGene) / (expected + rGene)`.
    #
    # rGene acts as a pseudo-count. When it is small relative to the gene's total
    # count, eta is determined by the data; when it is large, eta stays close to 1.
    "rGene": 20,


    # Prior strength for theta, the per-cell scale factor. theta is the ratio of
    # the cell's observed count to the count expected under a candidate cell type,
    # with a Gamma(rTheta, rTheta) prior of mean 1. The estimate is the posterior
    # mode, `theta = (observed + rTheta - 1) / (expected + rTheta)`.
    #
    # When rTheta is small relative to the cell's total count, theta is determined
    # by the data; when it is large, theta stays close to 1. Must be greater than
    # 1, otherwise a cell with no spots gets a theta of zero or less.
    #
    # A starting value is the average total count per cell, one number across all
    # classes. It is an indication only; the value that suits a dataset depends on
    # the data and can be much lower.
    #
    # A low value lets a cell with very few spots be claimed by a real class on
    # one or two spots; see
    # [the prior strength](../the-model/scale-factors.md#r-theta)
    # for what that costs.
    "rTheta": 25.0,


    # Dispersion of the negative binomial model for gene counts. The negative
    # binomial is used instead of a Poisson because gene counts vary more between
    # cells of the same type than a Poisson allows. rSpot is also the prior
    # strength for gamma, the per-cell, per-gene scale factor, which has a
    # Gamma(rSpot, rSpot) prior; integrating gamma out gives the negative binomial.
    # Lower values allow more variation in a gene's count between cells of the
    # same type; higher values allow less. The model approaches a Poisson only for
    # genes whose expected count in a cell is well below rSpot; for the genes a
    # class expects most, that would take a value in the hundreds.
    "rSpot": 2,


    # Bonus added to the log score of the cell whose segmentation label contains
    # the spot. 0 disables it. True is converted to 2, with a warning. Non-zero
    # values make spot assignment depend on the exact position of the
    # segmentation boundaries.
    "InsideCellBonus": 0,


    # Strength of the spatial prior, a Markov random field on the cell types. When
    # a cell is scored, each cell type receives a bonus of mrf_beta times the
    # distance-weighted sum of the neighbouring cells' probabilities for that type.
    # 0 disables the spatial prior.
    #
    # The prior is most informative where cell types are spatially organised, for
    # example in cortical layers, and least informative where they are intermixed.
    "mrf_beta": 1.0,


    # Groups of cell types that the spatial prior does not distinguish. A
    # neighbour's probability for any type in a group counts towards every type in
    # the group, so all types in a group receive the same spatial bonus and are
    # separated by their gene counts alone. This is useful when a rare type occurs
    # inside a region dominated by a closely related type.
    #
    # A list of groups, each a list of two or more cell type names matching the
    # columns of `scRNAseq`, for example `[["037 DG Glut", "038 DG-PIR Ex IMN"]]`.
    # A cell type can belong to one group only. None means no groups.
    #
    # The spatial bonuses of the types in a group add up, so a large group can take
    # cells from types outside it. Groups should be kept small. See
    # [Pooling sister classes](../the-model/cell-class.md#pooling-sister-classes).
    "mrf_pooled_classes": None,


    # Prior mean of the background density, in spots per unit volume. Spots not
    # explained by any cell, such as RNA in cell processes or technical misreads,
    # are modelled as a uniform background over the imaged region, and each spot
    # is assigned to one of its neighbouring cells or to the background. Higher
    # values assign more spots to the background; lower values assign more spots
    # to cells.
    #
    # A separate density is estimated for each gene, with a
    # Gamma(rRho, rRho / MisreadDensity) prior.
    #
    # Volume is measured in xy pixels, with z scaled by voxel_size[2] /
    # voxel_size[0]. For example, 100 x 100 pixels over 10 planes with voxel size
    # [0.28, 0.28, 0.7] is a volume of 100 x 100 x 10 x 2.5 = 250,000, and one
    # expected misread of a gene in that block is a density of 4e-6.
    #
    # A number or a dictionary with a 'default' key, for example
    # `{'default': 1e-6}`. A number is converted to `{'default': value}`. Only the
    # 'default' key is used; per-gene entries are ignored.
    "MisreadDensity": 0.00001,


    # Prior strength for the per-gene background density. The prior is
    # Gamma(rRho, rRho / MisreadDensity), with mean MisreadDensity, and rRho acts
    # as a pseudo-count against the number of spots of the gene assigned to the
    # background. When rRho is small relative to that number, the density is
    # determined by the data; when it is large, every gene stays close to
    # MisreadDensity.
    "rRho": 1.0,


    # Pseudo-count of the segmentation centroid in the cell geometry update. The
    # segmentation centroid counts for this many spots when it is averaged with the
    # spots assigned to the cell to give the fitted cell centroid.
    "cell_centroid_prior": 1000,


    # Degrees of freedom of the inverse Wishart prior on each cell's covariance.
    # The prior shape, a ball of the mean cell radius, counts for
    # cell_cov_prior - d - 1 spots, d being 2 or 3. Must be larger than d + 1.
    "cell_cov_prior": 1000,


    # Constant added to the expected counts, in spots per cell. Without it, one
    # spot of a gene with zero expression in a cell type would give that type zero
    # likelihood, however well the other genes match. The constant allows for such
    # spots arising from technical errors. It is also the expected count of every
    # gene in the Zero class.
    "SpotReg": 0.1,


    # Number of nearest cells considered for each spot, and number of neighbouring
    # cells used by the spatial prior. Each spot is scored against its nNeighbors
    # nearest cells and the background, so at 9 there are ten candidates. Higher
    # values let a spot be assigned to a more distant cell, at a higher cost per
    # iteration. Lower values are faster, but a spot near a cell boundary may miss
    # the cell it belongs to.
    "nNeighbors": 9,


    # Write the results to `<output_path>/pciSeq/data/`: cell, spot and cell
    # boundary tables as tsv and feather files, a SpatialData zarr store, the
    # diagnostics database, and the fitted model as a pickle.
    "save_data": True,


    # Log the run time of each step and, for each iteration, the cells and spots
    # with the largest changes. Does not affect the result, but adds computation
    # to every iteration.
    "verbose": False,


    # Steps of the variational loop at which the ELBO is evaluated, before and
    # after the step, with the change logged. Valid names are geneCount_upd,
    # rho_upd, eta_upd, theta_upd, gamma_upd, cell_to_cellType, dalpha_upd and
    # spots_to_cell; "all" selects every step. An empty list disables it.
    # Evaluating the ELBO is expensive, so selecting many steps slows the run
    # considerably. Intended for development.
    "elbo_per_step": [],


    # Directory for the results. Output is written to `<output_path>/pciSeq/`.
    # 'default' uses the system temporary directory.
    "output_path": "default",


    # Radius of the prior on each cell's Gaussian, in the units of the spot
    # coordinates. The prior mean covariance is cell_radius^2 times the identity;
    # the fitted covariance is updated from the assigned spots (see cell_cov_prior).
    # None uses half the mean equivalent radius, sqrt(area / pi), of the segmented
    # cells.
    "cell_radius": None,


    # How the prior over cell types is computed:
    #
    # - 'uniform': cell_type_weights are used unchanged throughout.
    # - 'weighted': the Zero weight from cell_type_weights is fixed, and the
    #   weights of the other cell types are re-estimated at each iteration from
    #   the current cell type probabilities, through a Dirichlet posterior.
    "cell_type_prior": "uniform",


    # Prior probabilities of the cell types. Types not listed share the remaining
    # probability equally, so {"Zero": 0.5} gives Zero 0.5 and divides the other
    # 0.5 equally among the remaining types. None gives every type the same
    # weight. Names that do not match a cell type are ignored with a warning, and
    # values summing to more than 1 are normalised with a warning. There is no
    # 'default' key.
    #
    # The Zero class has zero expected expression for every gene. It takes cells
    # with zero or few reads that no real class explains; such a cell can be a real
    # cell of a type the gene panel does not target.
    "cell_type_weights": {"Zero": 0.5},



    # *******************************************************************************
    # 3D settings
    # *******************************************************************************

    # Physical size of a voxel as [x, y, z], all in the same unit. Coordinates are
    # rescaled by the ratios to the x size, so that distances are correct when the
    # z step differs from the pixel size. For example, [0.147, 0.147, 0.9] for
    # 0.147 um pixels and a 0.9 um z step. [1, 1, 1] is isotropic and applies to 2D
    # data.
    "voxel_size": [1, 1, 1],


    # Remove cells that appear on a single z plane, which are usually segmentation
    # artefacts. 3D data only.
    "remove_flat_cells": True,



    # *******************************************************************************
    # Live viewer, optional
    # *******************************************************************************

    # Show the cell typing in a browser while it runs. See
    # [Live viewer](./live-viewer.md).
    "realtime_viewer": False,


    # Port of the live viewer server.
    "realtime_viewer_port": 5001,


    # Radius at which the live viewer draws every cell. None draws each cell at its
    # own size.
    "realtime_viewer_fixed_radius": None,
}


# Set by pciSeq during a run and stored in the same dictionary as the settings
# above.
#
# - `is3D`: whether the segmentation has more than one plane.
# - `img_dim`: width, height and number of planes of the label image.
# - `label_map`: mapping used to renumber non-sequential cell labels, or None.
RUNTIME_KEYS = ('is3D', 'img_dim', 'label_map')

