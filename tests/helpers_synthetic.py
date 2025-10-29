import numpy as np
import pandas as pd


def generate_spheres_grid(
    radius: int, simulated_counts: dict, spacing_factor: int, rng: np.random.Generator
):
    """
    Place one sphere per class on a grid. Returns:
    - spheres: list of (z, y, x, radius, label_id, class_name)
    - shape: (Z, Y, X)
    """
    class_names = np.array(list(simulated_counts.keys()))
    n_cells = len(class_names)

    # Grid tiling
    grid_y = int(np.sqrt(n_cells)) or 1
    grid_x = int(np.ceil(n_cells / grid_y))
    spacing = spacing_factor * radius

    y_coords = np.arange(grid_y) * spacing + spacing // 2
    x_coords = np.arange(grid_x) * spacing + spacing // 2
    X, Y = np.meshgrid(x_coords, y_coords)
    X = X.flatten()[:n_cells]
    Y = Y.flatten()[:n_cells]

    # Random z per sphere
    z_coords = rng.integers(2 * radius, 4 * radius + 1, size=n_cells)

    labels = np.arange(1, n_cells + 1)
    spheres = np.column_stack(
        (z_coords, Y, X, np.full(n_cells, radius), labels, class_names[:, None])
    )

    max_y = grid_y * spacing
    max_x = grid_x * spacing
    max_z = 6 * radius
    shape = (max_z, max_y, max_x)

    spheres_list = [tuple(row) for row in spheres]
    return spheres_list, shape


def create_3d_label_image(shape, spheres):
    import numpy as np

    label_img = np.zeros(shape, dtype=np.uint16)
    for zc, yc, xc, r, label_id, _class_name in spheres:
        zc = np.int32(zc)
        yc = np.int32(yc)
        xc = np.int32(xc)
        r = np.int32(r)
        label_id = np.int32(label_id)
        z_min = max(zc - r, 0)
        z_max = min(zc + r + 1, shape[0])
        y_min = max(yc - r, 0)
        y_max = min(yc + r + 1, shape[1])
        x_min = max(xc - r, 0)
        x_max = min(xc + r + 1, shape[2])
        zz, yy, xx = np.ogrid[z_min:z_max, y_min:y_max, x_min:x_max]
        mask = ((zz - zc) ** 2 + (yy - yc) ** 2 + (xx - xc) ** 2) <= r**2
        sub = label_img[z_min:z_max, y_min:y_max, x_min:x_max]
        sub[mask] = label_id
        label_img[z_min:z_max, y_min:y_max, x_min:x_max] = sub
    return label_img


def sample_point_clouds(spheres, cell_gene_counts, rng: np.random.Generator):
    """Sample Gaussian points for each gene occurrence around each sphere center."""
    all_points = []
    class_names = list(cell_gene_counts.keys())
    counts_list = list(cell_gene_counts.values())
    for i, (zc, yc, xc, r, label_id, _class_name) in enumerate(spheres):
        zc = np.int32(zc)
        yc = np.int32(yc)
        xc = np.int32(xc)
        r = np.int32(r)
        label_id = np.int32(label_id)
        cls = class_names[i]
        counts = counts_list[i]
        genes = list(counts.keys())
        gene_counts = np.asarray(list(counts.values()), dtype=np.int64)
        if gene_counts.sum() == 0:
            continue
        collapsed = np.repeat(genes, gene_counts)
        n_rows = collapsed.shape[0]
        points = rng.normal(loc=[zc, yc, xc], scale=r, size=(n_rows, 3))
        labels = np.full((n_rows, 1), label_id, dtype=np.uint32)
        actual_class = np.full((n_rows, 1), cls)
        arr = np.column_stack((collapsed, points, labels, actual_class))
        all_points.append(arr)
    if not all_points:
        return pd.DataFrame(
            columns=["gene_name", "z", "y", "x", "label", "actual_class"]
        )
    all_points = np.vstack(all_points)
    df = pd.DataFrame(
        all_points, columns=["gene_name", "z", "y", "x", "label", "actual_class"]
    )
    df[["x", "y", "z"]] = df[["x", "y", "z"]].astype(np.float32)
    df["label"] = df["label"].astype(np.uint32)
    return df


def simulate_nb_matrices(
    expression_matrix: pd.DataFrame,
    rng: np.random.Generator,
    r: int = 2,
    inefficiency: float = 1.0,
    rGene: int = 20,
    num_simulations: int = 1,
):
    """Negative binomial simulation matching the approach in simulations repo."""
    mu_values = expression_matrix.values.astype(float)
    p_values = r / (r + mu_values)
    simulated_matrices = []
    nG = mu_values.shape[0]
    for _ in range(num_simulations):
        sim_counts = rng.negative_binomial(r, p_values)
        if inefficiency != 1:
            eta = rng.gamma(rGene, 1 / rGene, nG)
            sim_counts = inefficiency * eta[:, None] * sim_counts
        sim_df = pd.DataFrame(
            sim_counts, index=expression_matrix.index, columns=expression_matrix.columns
        )
        simulated_matrices.append(sim_df)
    return simulated_matrices
