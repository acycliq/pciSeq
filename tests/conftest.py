import os
from pathlib import Path
import numpy as np
import pandas as pd
import pytest


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
        "Inefficiency": 1.0,
        "SpotReg": 0.1,
        "rSpot": 2,
        "InsideCellBonus": 0,
        "MisreadDensity": 0.0,
        "nNeighbors": 6,
        "CellCallTolerance": 0.5,
        "voxel_size": [1, 1, 1],
    }
