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


# ------------------------------------------------------------------ #
# 7. Confusion matrix
# ------------------------------------------------------------------ #
def confusion_matrix(cellData, actual_labels, class_names=None):
    """Probability-weighted (n_classes x n_classes) confusion matrix.

    For each cell with true class k, adds its full posterior vector
    to row k. So diagonal entry [k,k] is the total posterior mass
    pciSeq put on the correct class for cells of class k.
    """
    if class_names is None:
        predicted_classes = set()
        for cls_list in cellData["ClassName"].values:
            predicted_classes.update(cls_list)
        class_names = sorted(set(actual_labels.values()) | predicted_classes)

    class_to_idx = {c: i for i, c in enumerate(class_names)}
    nK = len(class_names)
    cm = np.zeros((nK, nK), dtype=float)

    for cell_num, cls_list, prob_list in zip(
        cellData["Cell_Num"].astype(int).values,
        cellData["ClassName"].values,
        cellData["Prob"].values,
    ):
        truth = actual_labels.get(int(cell_num))
        if truth is None:
            continue
        truth_idx = class_to_idx.get(truth)
        if truth_idx is None:
            continue
        for cls_name, prob in zip(cls_list, prob_list):
            pred_idx = class_to_idx.get(cls_name)
            if pred_idx is not None:
                cm[truth_idx, pred_idx] += float(prob)

    out = pd.DataFrame(cm, index=class_names, columns=class_names)
    out.index.name = "truth"
    out.columns.name = "predicted"
    return out


def diagonal_accuracy(cm):
    """Average fraction of posterior mass on the correct class."""
    mat = cm.values if hasattr(cm, 'values') else cm
    row_sums = mat.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return float(np.trace(mat / row_sums) / mat.shape[0])


def plot_confusion_matrix(cm, title="Confusion matrix", normalize=True):
    """Interactive plotly heatmap of the confusion matrix."""
    import plotly.express as px

    mat = cm.copy()
    if normalize:
        row_sum = mat.sum(axis=1).replace(0, 1.0)
        mat = mat.div(row_sum, axis=0)

    fig = px.imshow(
        mat.values, x=list(mat.columns), y=list(mat.index),
        labels=dict(x="Predicted class", y="True class", color="Mass"),
        color_continuous_scale="Viridis", zmin=0.0, zmax=1.0 if normalize else None,
        aspect="equal",
    )
    fig.update_traces(
        xgap=0.25, ygap=0.25,
        hovertemplate="Truth: %{y}<br>Predicted: %{x}<br>Mass: %{z:.3f}<extra></extra>",
        colorbar=dict(len=0.5, y=0.5),
    )
    fig.update_layout(
        title=title, xaxis_tickangle=45,
        xaxis=dict(showgrid=True, gridwidth=1, gridcolor="lightgray"),
        yaxis=dict(showgrid=True, gridwidth=1, gridcolor="lightgray"),
        autosize=False, width=1000, height=1000,
    )
    return fig


# ------------------------------------------------------------------ #
# 6. Run pciSeq
# ------------------------------------------------------------------ #
def run_pciseq(spots_df, label_image, reference, rSpot=2.0, opts=None):
    """Feed spots and the label image into pciSeq.fit.

    Three pieces of plumbing:
      1. Rename z -> z_plane, add intensity and score columns
         (pciSeq expects them, our simulator doesn't have them).
      2. Slice the 3D label image into per-plane sparse matrices.
      3. Build the opts dict with defaults matching the simulator.

    rSpot is its own arg so the caller can keep it locked to the
    same value used in simulate_nb_cells.
    """
    from scipy.sparse import coo_matrix
    import pciSeq

    spots = spots_df.copy()
    spots = spots.rename(columns={"z": "z_plane"})
    spots["intensity"] = np.float32(1.0)
    spots["score"] = np.float32(1.0)

    coo = [coo_matrix(plane) for plane in label_image]

    if opts and "rSpot" in opts:
        raise ValueError("Pass rSpot as its own arg, not inside opts. It must stay locked to the simulator's value.")

    final_opts = {
        "Inefficiency": 1.0,
        "SpotReg": 0.1,
        "rSpot": rSpot,
        "InsideCellBonus": 0.0,
        "MisreadDensity": 3e-20,
        "nNeighbors": 6,
        "CellCallTolerance": 0.02,
        "voxel_size": [1, 1, 1],
        "rTheta": 2,  # mean total gene count per cell across all classes in the Yao reference
        "mrf_beta": 0,
    }
    if opts:
        final_opts.update(opts)

    pciSeq.attach_to_log()
    cellData, geneData = pciSeq.fit(spots=spots, coo=coo, scRNAseq=reference, opts=final_opts)
    return cellData, geneData


if __name__ == "__main__":
    import pciSeq
    pciSeq.attach_to_log()

    # ---- All settings in one place -------------------------------- #
    N_PER_CLASS    = 10000
    N_RUNS         = 100
    RSPOT          = 2.0
    RADIUS         = 18
    SPACING_FACTOR = 2
    SEED           = 42
    SHOW_NAPARI    = False
    # --------------------------------------------------------------- #

    # 1. get the single cell data
    scRNAseq = get_scRNAseq()
    logger.info(f"reference: {scRNAseq.shape[0]} genes x {scRNAseq.shape[1]} classes")

    # 2. simulate single cell data
    rng = np.random.default_rng(SEED)
    sim_dfs = simulate_nb_cells(scRNAseq, n_per_class=N_PER_CLASS, rSpot=RSPOT, rng=rng)
    logger.info(f"generated {len(sim_dfs)} sets, each {sim_dfs[0].shape}")

    # 3-7. loop: place on grid, make pointclouds, build label image, run pciSeq, accumulate confusion matrix
    GLOBAL_CLASSES = sorted(scRNAseq.columns.tolist())
    cm_total = pd.DataFrame(0.0, index=GLOBAL_CLASSES, columns=GLOBAL_CLASSES)
    cm_total.index.name = "truth"
    cm_total.columns.name = "predicted"

    for i in range(N_RUNS):
        logger.info(f"[run {i + 1}/{N_RUNS}]")
        cells_df = sim_dfs[i]

        # Shuffle which class lands at which grid position
        shuffled_cols = list(cells_df.columns)
        rng.shuffle(shuffled_cols)
        cells_df = cells_df[shuffled_cols]

        cell_grid = make_grid(cells_df, radius=RADIUS, spacing_factor=SPACING_FACTOR, rng=rng)
        spots_df = make_pointclouds(cells_df, cell_grid, radius=RADIUS, rng=rng)
        label_image = build_label_image(cell_grid, radius=RADIUS)
        cellData, geneData = run_pciseq(spots_df, label_image, scRNAseq, rSpot=RSPOT)
        actual_labels = dict(zip(cell_grid["cell_label"].astype(int), cell_grid["class_name"]))
        cm = confusion_matrix(cellData, actual_labels, class_names=GLOBAL_CLASSES)
        cm_total = cm_total + cm

    cm_avg = cm_total / N_RUNS
    diag_acc = diagonal_accuracy(cm_avg)
    logger.info(f"diagonal accuracy: {diag_acc:.4f} (1.0 = perfect, {1/len(GLOBAL_CLASSES):.3f} = chance)")

    fig = plot_confusion_matrix(cm_avg, title=f"Confusion matrix (diag acc = {diag_acc:.3f}, N={N_RUNS})")

    results_dir = Path(__file__).resolve().parent / "results" / "harness"
    results_dir.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(results_dir / "confusion_matrix.html"))
    logger.info(f"saved confusion matrix to {results_dir / 'confusion_matrix.html'}")
    fig.show()

    # View in napari (uses the last iteration's data)
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