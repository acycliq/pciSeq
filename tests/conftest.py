import os
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from pciSeq.src.core.main import VarBayes


@pytest.fixture(scope="session")
def rng():
    return np.random.default_rng(42)


@pytest.fixture(scope="session")
def scref_df():
    # Prefer vendored data, fall back to env var
    vendored = Path(__file__).parent / "data" / "scRNAseq_final.csv"
    env_path = os.environ.get("SCREF_PATH")
    if vendored.exists():
        df = pd.read_csv(vendored)
    elif env_path and Path(env_path).exists():
        df = pd.read_csv(env_path)
    else:
        pytest.skip(
            "Missing scRNA reference CSV (tests/data/scRNAseq_final.csv or SCREF_PATH)"
        )

    # Expect index column named 'Unnamed: 0' like in reference_data.py
    if "Unnamed: 0" in df.columns:
        df = df.set_index("Unnamed: 0")
    # Set axis names for clarity
    df.index.name = "gene_name"
    df.columns.name = "class_name"
    # Remove duplicate genes by keeping max total, like reference_data.keep_labels_unique
    if df.index.duplicated().any():
        tmp = df.copy().assign(total=df.sum(axis=1))
        tmp = (
            tmp.sort_values(["gene_name", "total"], ascending=[True, False])
            .groupby("gene_name")
            .head(1)
            .drop(columns=["total"])
        )
        df = tmp
    # Drop all-zero columns
    df = df.loc[:, (df != 0).any(axis=0)]
    return df


@pytest.fixture(scope="session")
def base_opts():
    # Match main.py defaults for Case 1 (no inefficiency, no misreads)
    return {
        "exclude_genes": [],
        "max_iter": 1000,
        "CellCallTolerance": 0.5,
        "rGene": 20,
        "Inefficiency": 1.0,
        "InsideCellBonus": 0,
        "mrf_beta": 1.0,
        "MisreadDensity": {"default": 1e-6},
        "cell_centroid_prior": 10,
        "cell_cov_prior": 10,
        "SpotReg": 0.1,
        "nNeighbors": 6,
        "rSpot": 2,
        "save_data": False,
        "output_path": "default",
        "launch_viewer": False,
        "launch_diagnostics": False,
        "is_redis_running": False,
        "cell_radius": None,
        "cell_type_prior": "uniform",
        "cell_type_weights": None,
        "is3D": True,
        "voxel_size": [1, 1, 1],
        "exclude_planes": [],
        "remove_flat_cells": True,
        "mean_gene_counts_per_class": 60,
        "mean_gene_counts_per_cell": 30,
        "img_dim": {"w": 100, "h": 100, "n_planes": 10},
        "rRho": 1000.0,
        "rTheta": 25.0,
        "label_map": {},
    }


@pytest.fixture
def minimal_varbayes(rng, base_opts):
    """
    Create a minimal VarBayes instance for unit testing.

    Small synthetic data:
    - 5 cells (including background cell 0)
    - 10 genes
    - 3 cell types
    - 30 spots
    """
    # Create minimal scRNA reference (10 genes x 3 cell types)
    gene_names = [f"Gene_{i}" for i in range(10)]
    cell_types = ["Type_A", "Type_B", "Type_C"]

    # Create expression profiles with some structure
    scref_data = rng.random((10, 3)) * 100
    scref_df = pd.DataFrame(scref_data, index=gene_names, columns=cell_types)
    scref_df.index.name = "gene_name"
    scref_df.columns.name = "class_name"

    # Create minimal cells DataFrame (8 cells - need more than nNeighbors=6)
    # Columns match what stage_data produces: label, area, x0, y0, z0, values
    cells_df = pd.DataFrame(
        {
            "label": [
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
            ],  # No cell 0 - background is handled separately
            "area": [100, 120, 90, 110, 95, 105, 115, 100],  # Cell areas
            "x0": [10, 20, 30, 40, 10, 20, 30, 40],  # Cell centroids
            "y0": [10, 20, 30, 40, 50, 50, 50, 50],
            "z0": [5, 5, 5, 5, 5, 5, 5, 5],
            "values": [
                [[10, 10, 5]],
                [[20, 20, 5]],
                [[30, 30, 5]],
                [[40, 40, 5]],
                [[10, 50, 5]],
                [[20, 50, 5]],
                [[30, 50, 5]],
                [[40, 50, 5]],
            ],  # Cell voxel coordinates
        }
    )

    # Create minimal spots DataFrame (30 spots)
    # Columns match what stage_data produces: x, y, z, plane_id, label, gene_name, score, intensity
    spots_df = pd.DataFrame(
        {
            "x": rng.uniform(0, 50, 30),
            "y": rng.uniform(0, 50, 30),
            "z": rng.uniform(0, 10, 30),
            "plane_id": rng.integers(0, 3, 30),
            "label": rng.choice([1, 2, 3, 4, 5, 6, 7, 8], 30),  # Assign spots to cells
            "gene_name": rng.choice(gene_names, 30),
            "score": rng.uniform(0.5, 1.0, 30),
            "intensity": rng.uniform(100, 1000, 30),
        }
    )

    # Instantiate VarBayes
    vb = VarBayes(
        spots_df=spots_df, cells_df=cells_df, scRNAseq=scref_df, config=base_opts
    )

    return vb
