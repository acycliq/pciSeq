import numpy as np
import pytest
from scipy.sparse import coo_matrix

import pciSeq
from .helpers_synthetic import (
    generate_spheres_grid,
    create_3d_label_image,
    sample_point_clouds,
    simulate_nb_matrices,
)


@pytest.mark.slow
def test_case1_spacing4_full_panel_high_accuracy(rng, scref_df, base_opts):
    radius = 18
    spacing_factor = 4

    # One simulation of full panel counts with NB r = rSpot
    sim_df = simulate_nb_matrices(
        scref_df,
        rng,
        r=base_opts["rSpot"],
        inefficiency=base_opts["Inefficiency"],
        num_simulations=1,
    )[0]

    # Build dict class -> {gene: count}
    simulated_counts = sim_df.to_dict()

    # Shuffle class order to avoid positional bias
    keys = list(simulated_counts.keys())
    rng.shuffle(keys)
    simulated_counts = {k: simulated_counts[k] for k in keys}

    # Geometry and labels
    spheres, shape = generate_spheres_grid(
        radius, simulated_counts, spacing_factor, rng
    )
    label_img = create_3d_label_image(shape, spheres)

    # Spots
    spots = sample_point_clouds(spheres, simulated_counts, rng)
    if spots.empty:
        pytest.skip("No spots generated (all-zero simulation)")

    # add summy values for score and intensity
    spots = spots.assign(score=0)
    spots = spots.assign(intensity=0)
    spots = spots.rename(columns={"z": "z_plane", "gene_name": "gene_name"})

    # coo per z-slice
    coo_slices = [coo_matrix(plane) for plane in label_img]

    # Run pciSeq
    pciSeq.attach_to_log()
    cellData, geneData = pciSeq.fit(
        spots=spots, coo=coo_slices, scRNAseq=scref_df, opts=base_opts
    )

    # Basic structure checks
    assert cellData is not None and geneData is not None
    assert {"Cell_Num", "ClassName", "Prob"}.issubset(
        set(cellData.columns)
    ), "Missing columns in cellData"
    assert cellData.shape[0] == len(
        spheres
    ), "cellData rows must equal number of spheres"

    # Label alignment: raster labels are 1..N, check presence
    raster_labels = sorted([int(s[4]) for s in spheres])
    assert sorted(cellData.Cell_Num.tolist()) == raster_labels

    # Compute accuracy
    cell_label = cellData.Cell_Num.values
    best_class = cellData.ClassName.map(lambda d: d[0]).values
    actual_class = np.array([s[5] for s in spheres])

    # Sort by Cell_Num to align
    order = np.argsort(cell_label)
    best_class = best_class[order]
    actual_class = actual_class[np.argsort([int(s[4]) for s in spheres])]

    accuracy = float(np.mean(best_class == actual_class))

    # Expect high accuracy for spacing_factor=4 without inefficiency or misreads
    assert (
        accuracy >= 0.9
    ), f"Accuracy {accuracy:.3f} below expected threshold for spacing=4"
