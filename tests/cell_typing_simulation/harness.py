"""
Simulate gene counts from a Negative Binomial, place them as point clouds
inside spheres on a 3D grid, feed everything to pciSeq and check whether
pciSeq correctly recovers the class that generated the counts.

Built step-by-step, following the same logic as the original harness at
/home/dimitris/dev/python/negative_binomial_simulations/main.py
"""
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

logger = logging.getLogger(__name__)


DATA_DIR = Path(__file__).resolve().parent / "data" / "yao"


# ------------------------------------------------------------------ #
# 1. Load the single cell reference data
# ------------------------------------------------------------------ #
def get_scRNAseq(path=None):
    if path is None:
        path = DATA_DIR / "scRNAseq_final.csv"
    scRNAseq = pd.read_csv(path).set_index("Unnamed: 0")
    scRNAseq = scRNAseq.rename_axis("class_name", axis="columns").rename_axis("gene_name")

    # Drop duplicate gene labels, keep the one with highest total
    scRNAseq = scRNAseq.assign(total=scRNAseq.sum(axis=1))
    scRNAseq = (scRNAseq.sort_values(["gene_name", "total"], ascending=[True, False])
                        .groupby("gene_name").head(1))
    scRNAseq = scRNAseq.drop(columns=["total"])

    # Drop classes that are zero everywhere
    scRNAseq = scRNAseq.loc[:, (scRNAseq != 0).any(axis=0)]

    # Aggregate sub-classes by name
    scT = scRNAseq.T
    scRNAseq = scT.groupby(scT.index.values).agg("mean").T
    return scRNAseq


# ------------------------------------------------------------------ #
# 2. Simulate gene counts from the Negative Binomial
# ------------------------------------------------------------------ #
def simulate_nb_cells(reference, n_per_class=1000, rSpot=2.0, rng=None):
    """For each class, draw n_per_class independent NB cells.

    Returns a list of DataFrames (one per simulation), each with the
    same shape as the reference: (n_genes x n_classes).
    """
    if rng is None:
        rng = np.random.default_rng()

    mu = reference.values.astype(float)
    p = rSpot / (rSpot + mu)

    # One big vectorised draw: (n_per_class, n_genes, n_classes)
    counts = rng.negative_binomial(rSpot, p, size=(n_per_class, *p.shape))

    return [
        pd.DataFrame(counts[i], index=reference.index, columns=reference.columns)
        for i in range(n_per_class)
    ]


# ------------------------------------------------------------------ #
# 3. Place cells on a grid
# ------------------------------------------------------------------ #
def make_grid(cells_df, radius=18, spacing_factor=6, rng=None):
    """Place one sphere per class on a regular xy grid with random z.

    Takes one DataFrame from simulate_nb_cells (n_genes x n_classes).
    Each column becomes one sphere. Returns a DataFrame with one row
    per placed cell: cell_label, class_name, z, y, x.
    """
    if rng is None:
        rng = np.random.default_rng()

    class_names = list(cells_df.columns)
    n_cells = len(class_names)

    grid_y = int(np.sqrt(n_cells))
    grid_x = int(np.ceil(n_cells / grid_y))
    spacing = spacing_factor * radius

    y_coords = np.arange(grid_y) * spacing + spacing // 2
    x_coords = np.arange(grid_x) * spacing + spacing // 2
    X, Y = np.meshgrid(x_coords, y_coords)
    X = X.flatten()[:n_cells]
    Y = Y.flatten()[:n_cells]

    z_coords = rng.integers(2 * radius, 4 * radius + 1, size=n_cells)

    return pd.DataFrame({
        "cell_label": np.arange(n_cells, dtype=np.int32) + 1,
        "class_name": class_names,
        "z": z_coords.astype(np.int32),
        "y": Y.astype(np.int32),
        "x": X.astype(np.int32),
    })


# ------------------------------------------------------------------ #
# 4. Sample point clouds inside each sphere
# ------------------------------------------------------------------ #
def make_pointclouds(cells_df, cell_grid, radius=18, rng=None):
    """For each cell on the grid, draw spots from a 3D Gaussian.

    For each placed cell, looks up its gene counts in cells_df and
    draws that many spots from N(centroid, radius * I). Returns a
    long-form DataFrame: gene_name, z, y, x, cell_label, class_name.
    """
    if rng is None:
        rng = np.random.default_rng()

    pointclouds = []
    for row in cell_grid.itertuples(index=False):
        counts_vec = cells_df[row.class_name].values
        gene_names = cells_df.index.values
        n_spots = int(counts_vec.sum())
        if n_spots == 0:
            continue

        positions = rng.normal(
            loc=[row.z, row.y, row.x],
            scale=radius,
            size=(n_spots, 3),
        )
        pointclouds.append(pd.DataFrame({
            "gene_name": np.repeat(gene_names, counts_vec),
            "z": positions[:, 0].astype(np.float32),
            "y": positions[:, 1].astype(np.float32),
            "x": positions[:, 2].astype(np.float32),
            "cell_label": np.int32(row.cell_label),
            "class_name": row.class_name,
        }))

    if not pointclouds:
        return pd.DataFrame(columns=["gene_name", "z", "y", "x",
                                      "cell_label", "class_name"])
    return pd.concat(pointclouds, ignore_index=True)


# ------------------------------------------------------------------ #
# 5. Build the 3D label image
# ------------------------------------------------------------------ #
def build_label_image(cell_grid, radius=18):
    """Paint each cell as a solid sphere in a uint16 3D volume.

    Background voxels are 0, voxels inside cell k carry cell_grid.cell_label[k].
    Image bounds are derived from the grid with a 2*radius margin.
    """
    margin = 2 * radius
    n_z = int(cell_grid["z"].max() + margin)
    n_y = int(cell_grid["y"].max() + margin)
    n_x = int(cell_grid["x"].max() + margin)

    label_image = np.zeros((n_z, n_y, n_x), dtype=np.uint16)
    r = int(radius)

    for row in cell_grid.itertuples(index=False):
        zc, yc, xc = int(row.z), int(row.y), int(row.x)
        lbl = int(row.cell_label)

        z_min, z_max = max(zc - r, 0), min(zc + r + 1, n_z)
        y_min, y_max = max(yc - r, 0), min(yc + r + 1, n_y)
        x_min, x_max = max(xc - r, 0), min(xc + r + 1, n_x)

        zz, yy, xx = np.ogrid[z_min:z_max, y_min:y_max, x_min:x_max]
        mask = (zz - zc) ** 2 + (yy - yc) ** 2 + (xx - xc) ** 2 <= r ** 2

        sub = label_image[z_min:z_max, y_min:y_max, x_min:x_max]
        sub[mask] = lbl

    return label_image


if __name__ == "__main__":
    import pciSeq
    pciSeq.attach_to_log()

    N_PER_CLASS = 10000
    RSPOT = 2.0
    SEED = 42

    # 1. get the single cell data
    scRNAseq = get_scRNAseq()
    logger.info(f"reference: {scRNAseq.shape[0]} genes x {scRNAseq.shape[1]} classes")

    # 2. simulate single cell data
    rng = np.random.default_rng(SEED)
    sim_dfs = simulate_nb_cells(scRNAseq, n_per_class=N_PER_CLASS, rSpot=RSPOT, rng=rng)
    logger.info(f"generated {len(sim_dfs)} sets, each {sim_dfs[0].shape}")

    # Sanity check: the mean across all sets should be close to the reference
    avg = sum(sim_dfs) / len(sim_dfs)
    cls = scRNAseq.columns[0]
    top_genes = scRNAseq[cls].nlargest(5).index
    logger.info(f"class: {cls}")
    logger.info(f"  reference mean (top 5 genes): {scRNAseq[cls][top_genes].to_dict()}")
    logger.info(f"  simulated mean (top 5 genes): {avg[cls][top_genes].round(2).to_dict()}")

    # 3. place cells on a grid
    RADIUS = 18
    SPACING_FACTOR = 2
    cells_df = sim_dfs[0]
    cell_grid = make_grid(cells_df, radius=RADIUS, spacing_factor=SPACING_FACTOR, rng=rng)
    logger.info(f"placed {len(cell_grid)} cells on a grid "
                f"(radius={RADIUS}, spacing={SPACING_FACTOR}x)")
    logger.info(f"\n{cell_grid.head(5).to_string(index=False)}")

    # 4. sample point clouds inside each sphere
    spots_df = make_pointclouds(cells_df, cell_grid, radius=RADIUS, rng=rng)
    logger.info(f"sampled {len(spots_df)} spots across {cell_grid.cell_label.nunique()} cells")
    logger.info(f"  columns: {list(spots_df.columns)}")
    logger.info(f"  spots per cell (first 5): "
                f"{spots_df.groupby('cell_label').size().head(5).to_dict()}")

    # 5. build the 3D label image
    label_image = build_label_image(cell_grid, radius=RADIUS)
    logger.info(f"label image shape (z, y, x) = {label_image.shape}, "
                f"nonzero voxels = {int((label_image > 0).sum())}")

    # View in napari
    SHOW_NAPARI = True
    if SHOW_NAPARI:
        import napari
        viewer = napari.Viewer(ndisplay=3)
        viewer.add_labels(label_image, name="Cells")
        viewer.add_points(
            spots_df[["z", "y", "x"]].values,
            properties={
                "label": spots_df["cell_label"].values,
                "gene": spots_df["gene_name"].values,
            },
            face_color="label",
            face_colormap="turbo",
            size=1.5,
            name="Spots",
        )
        napari.run()