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
    SPACING_FACTOR = 6
    cells_df = sim_dfs[0]
    cell_grid = make_grid(cells_df, radius=RADIUS, spacing_factor=SPACING_FACTOR, rng=rng)
    logger.info(f"placed {len(cell_grid)} cells on a grid "
                f"(radius={RADIUS}, spacing={SPACING_FACTOR}x)")
    logger.info(f"\n{cell_grid.head(5).to_string(index=False)}")