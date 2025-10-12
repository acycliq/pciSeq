"""File and data I/O operations."""
import os
import numpy as np
from pathlib import Path
import pickle
import tempfile
import shutil
import pyarrow as pa
import pyarrow.feather as feather
from pathlib import Path
import json
from typing import Tuple, Optional, Dict, Any, Union
from urllib.parse import urlparse
from urllib.request import urlopen
from typing import List, Any, Dict
import pandas as pd
from tqdm import tqdm
import logging

# Configure logging
io_utils_logger = logging.getLogger(__name__)


def get_out_dir(path: Optional[str] = None, sub_folder: str = '') -> str:
    """Get or create output directory path.

    Args:
        path: Base path, or None for default temp directory
        sub_folder: Optional subdirectory name

    Returns:
        str: Path to output directory

    Notes:
        - Uses system temp directory if path is None or 'default'
        - Creates directories if they don't exist
    """
    if path is None or path == 'default':
        out_dir = Path(tempfile.gettempdir()) / 'pciSeq'
    else:
        out_dir = Path(path) / sub_folder / 'pciSeq'

    out_dir.mkdir(parents=True, exist_ok=True)
    return str(out_dir)


def log_file(cfg: Dict) -> None:
    """Setup the logger file handler if it doesn't exist.

    Args:
        cfg: Configuration dictionary containing output path

    Notes:
        - Only adds FileHandler if one doesn't already exist
        - Creates log file in output directory
        - Uses standard logging format
    """
    root_logger = logging.getLogger()

    # Only add FileHandler if none exists
    if root_logger.handlers and not any(isinstance(h, logging.FileHandler) for h in root_logger.handlers):
        logfile = os.path.join(get_out_dir(cfg['output_path']), 'pciSeq.log')
        fh = logging.FileHandler(logfile, mode='w')
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)

        root_logger.addHandler(fh)
        io_utils_logger.info('Writing to %s' % logfile)


def download_url_to_file(url: str, dst: str, progress: bool = True) -> None:
    """Download object at the given URL to a local path.

    Args:
        url: URL of the object to download
        dst: Full path where object will be saved
        progress: Whether to display a progress bar

    Notes:
        Thanks to torch, slightly modified
    """
    file_size = None
    u = urlopen(url)
    meta = u.info()

    if hasattr(meta, 'getheaders'):
        content_length = meta.getheaders("Content-Length")
    else:
        content_length = meta.get_all("Content-Length")

    if content_length is not None and len(content_length) > 0:
        file_size = int(content_length[0])

    # Save to temp file first then move
    dst = os.path.expanduser(dst)
    dst_dir = os.path.dirname(dst)
    f = tempfile.NamedTemporaryFile(delete=False, dir=dst_dir)

    try:
        with tqdm(total=file_size, disable=not progress,
                  unit='B', unit_scale=True, unit_divisor=1024) as pbar:
            while True:
                buffer = u.read(8192)
                if len(buffer) == 0:
                    break
                f.write(buffer)
                pbar.update(len(buffer))
        f.close()
        shutil.move(f.name, dst)
    finally:
        f.close()
        if os.path.exists(f.name):
            os.remove(f.name)


def load_from_url(url: str) -> str:
    """Download file from URL if not already present locally.

    Args:
        url: URL to download from

    Returns:
        str: Local filename
    """
    parts = urlparse(url)
    filename = os.path.basename(parts.path)
    if not os.path.exists(filename):
        io_utils_logger.info('Downloading: "%s" to %s', url, filename)
        download_url_to_file(url, filename)
    return filename


def serialise(varBayes: Any, debug_dir: str) -> None:
    """Pickle variable Bayes object to debug directory.

    Args:
        varBayes: Object to serialize
        debug_dir: Directory to save pickle file
    """
    if not os.path.exists(debug_dir):
        os.makedirs(debug_dir)
    pickle_dst = os.path.join(debug_dir, 'pciSeq.pickle')
    with open(pickle_dst, 'wb') as outf:
        pickle.dump(varBayes, outf)
        io_utils_logger.info('Saved at %s', pickle_dst)


def export_db_tables(out_dir: str, con: Any) -> None:
    """Export all database tables to CSV files.

    Args:
        out_dir: Output directory for CSV files
        con: Database connection object
    """
    tables = con.get_db_tables()
    for table in tables:
        export_db_table(table, out_dir, con)


def export_db_table(table_name: str, out_dir: str, con: Any) -> None:
    """Export single database table to CSV.

    Args:
        table_name: Name of table to export
        out_dir: Output directory
        con: Database connection object
    """
    df = con.from_redis(table_name)
    fname = os.path.join(out_dir, table_name + '.csv')
    df.to_csv(fname, index=False)
    io_utils_logger.info('Saved at %s', fname)


def write_data(cellData: pd.DataFrame, geneData: pd.DataFrame,
               cellBoundaries: pd.DataFrame, cellBoundaries_list: pd.DataFrame, varBayes: Any, cfg: Dict) -> None:

    dst = get_out_dir(cfg['output_path'])
    out_dir = os.path.join(dst, 'data')

    write_tsv(cellData, geneData, cellBoundaries, out_dir)
    write_arrow(geneData, cellData, cellBoundaries_list, out_dir)

    # Save debug info
    serialise(varBayes, os.path.join(out_dir, 'debug'))


def write_tsv(cellData: pd.DataFrame, geneData: pd.DataFrame, cellBoundaries: pd.DataFrame,
              out_dir: os.path) -> None:
    """Write all data files to output directory.

    Args:
        cellData: Cell data DataFrame
        geneData: Gene data DataFrame
        cellBoundaries: Cell boundaries DataFrame
        out_dir: path to output directory
    """
    out_dir = os.path.join(out_dir, 'tsv')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Save cell data
    cellData.to_csv(os.path.join(out_dir, 'cellData.tsv'), sep='\t', index=False)
    io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'cellData.tsv'))

    # Save gene data
    geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
    io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'geneData.tsv'))

    # Save boundaries
    cellBoundaries.to_csv(os.path.join(out_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
    io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'cellBoundaries.tsv'))


def write_arrow(geneData:pd.DataFrame, cellData:pd.DataFrame, cellBoundaries:pd.DataFrame, out_dir: str = None) -> None:
    geneData_to_arrow(geneData, out_dir)
    cellData_to_arrow(cellData, out_dir)
    # io_utils_logger.info('boundaries_to_arrow_old - Starting')
    # boundaries_to_arrow_old(cellBoundaries, out_dir)
    # io_utils_logger.info('boundaries_to_arrow_old - Ending')

    # io_utils_logger.info('boundaries_to_arrow - Starting')
    boundaries_to_arrow(cellBoundaries, out_dir)
    # io_utils_logger.info('boundaries_to_arrow - Ending')

    # io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'cellBoundaries.tsv'))


def geneData_to_arrow(df_in: pd.DataFrame, out_dir: str = None) -> None:

    out_dir = Path(out_dir) / "arrow" / 'arrow_spots'
    out_dir.mkdir(parents=True, exist_ok=True)

    shards = []
    total_rows = 0
    shard_index = 0
    gene_dict_data = {}

    chunk_size = 200000
    for start in range(0, len(df_in), chunk_size):
        df = df_in.iloc[start:start+chunk_size]

        # Build gene dict from this chunk
        chunk_gene_dict = dict(zip(df["gene_id"], df["gene_name"]))
        gene_dict_data.update(chunk_gene_dict)

        # Create arrays with exact schema matching working converter
        # Column order: ['x', 'y', 'z', 'plane_id', 'spot_id', 'gene_id', 'neighbour_array', 'neighbour_prob', 'omp_score', 'omp_intensity']
        arrays = {}

        # Required columns with exact types from working converter
        if "x" in df.columns:
            arrays["x"] = pa.array(df["x"].astype("float32"))
        if "y" in df.columns:
            arrays["y"] = pa.array(df["y"].astype("float32"))
        if "z" in df.columns:
            arrays["z"] = pa.array(df["z"].astype("float32"))
        if "plane_id" in df.columns:
            arrays["plane_id"] = pa.array(df["plane_id"].astype("uint16"))
        if "spot_id" in df.columns:
            arrays["spot_id"] = pa.array(df["spot_id"].astype("uint32"))
        if "gene_id" in df.columns:
            arrays["gene_id"] = pa.array(df["gene_id"].astype("uint32"))

        # List columns
        if "neighbour_array" in df.columns:
            arrays["neighbour_array"] = pa.array(df["neighbour_array"].tolist(), type=pa.list_(pa.int32()))
        if "neighbour_prob" in df.columns:
            arrays["neighbour_prob"] = pa.array(df["neighbour_prob"].tolist(), type=pa.list_(pa.float32()))

        # Optional OMP columns
        if "omp_score" in df.columns:
            arrays["omp_score"] = pa.array(df["omp_score"].astype("float32"))
        if "omp_intensity" in df.columns:
            arrays["omp_intensity"] = pa.array(df["omp_intensity"].astype("float32"))

        # NOTE: gene_name and neighbour columns are excluded to match working converter

        table = pa.table(arrays)
        shard_name = f"spots_shard_{shard_index:03d}.feather"
        shard_path = out_dir / shard_name
        feather.write_feather(table, shard_path.as_posix(), compression='uncompressed')

        row_count = len(df)
        shards.append({"url": shard_name, "rows": int(row_count)})
        total_rows += row_count
        shard_index += 1

    # Write manifest
    manifest = {
        "format": "arrow-feather",
        "total_rows": int(total_rows),
        "shards": shards,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    # Write gene dictionary (id -> name) using data collected during chunking
    (out_dir / "gene_dict.json").write_text(json.dumps(gene_dict_data, indent=2))

    # io_utils_logger.info(f"Saved {total_rows} rows in {len(shards)} shards at {out_dir}")
    io_utils_logger.info(f"Saved at {out_dir}")



def cellData_to_arrow(df_in: pd.DataFrame, out_dir: str = None) -> None:
    out_dir = Path(out_dir) / "arrow" / 'arrow_cells'
    out_dir.mkdir(parents=True, exist_ok=True)

    shards = []
    total_rows = 0
    shard_index = 0

    chunk_size = 100000  # Match the working cell converter chunk size
    for start in range(0, len(df_in), chunk_size):
        df = df_in.iloc[start:start+chunk_size]

        # Create arrays with exact schema matching working cell converter
        # Expected columns: ['cell_id', 'X', 'Y', 'Z', 'class_name', 'prob', 'gaussian_contour', 'sphere_scale', 'sphere_rotation']
        arrays = {}

        # Map Cell_Num to cell_id with int32 type
        if "Cell_Num" in df.columns:
            arrays["cell_id"] = pa.array(df["Cell_Num"].astype("int32"))

        # Coordinate columns as float32
        if "X" in df.columns:
            arrays["X"] = pa.array(df["X"].astype("float32"))
        if "Y" in df.columns:
            arrays["Y"] = pa.array(df["Y"].astype("float32"))
        if "Z" in df.columns:
            arrays["Z"] = pa.array(df["Z"].astype("float32"))

        # List columns for cell classification - map from ClassName/Prob to class_name/prob
        if "ClassName" in df.columns:
            arrays["class_name"] = pa.array(df["ClassName"].tolist(), type=pa.list_(pa.string()))
        if "Prob" in df.columns:
            arrays["prob"] = pa.array(df["Prob"].tolist(), type=pa.list_(pa.float32()))

        # String columns (not nested lists like the bad files had)
        if "gaussian_contour" in df.columns:
            arrays["gaussian_contour"] = pa.array(df["gaussian_contour"].astype("string"))
        if "sphere_scale" in df.columns:
            arrays["sphere_scale"] = pa.array(df["sphere_scale"].astype("string"))
        if "sphere_rotation" in df.columns:
            arrays["sphere_rotation"] = pa.array(df["sphere_rotation"].astype("string"))

        # NOTE: Exclude columns that aren't in the working schema:
        # - Genenames, CellGeneCount, spot_id (these are not in the working cell files)

        table = pa.table(arrays)
        shard_name = f"cells_shard_{shard_index:03d}.feather"
        feather.write_feather(table, (out_dir / shard_name).as_posix(), compression='uncompressed')
        n = len(df)
        shards.append({"url": shard_name, "rows": int(n)})
        total_rows += n
        shard_index += 1

    # Manifest only - no class dict needed since class names are in the Arrow files
    manifest = {"format": "arrow-feather", "total_rows": int(total_rows), "shards": shards}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    # io_utils_logger.info(f"Saved {total_rows} rows in {len(shards)} shards at {out_dir}")
    io_utils_logger.info(f"Saved at {out_dir}")


def parse_coords(cell: str) -> List[Tuple[float, float]]:
    """Parse a JSON-like coords string '[[x,y], ...]' into a list of (x,y)."""
    if pd.isna(cell):
        return []
    try:
        arr = json.loads(cell)
        # Expect list of [x, y]
        out: List[Tuple[float, float]] = []
        for pair in arr:
            if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                continue
            x, y = pair
            out.append((float(x), float(y)))
        return out
    except Exception:
        return []


def validate_df_structure(df):
    """Validate that dataframe has required columns: plane_id, label, coords"""
    required_columns = {"plane_id", "label", "coords"}

    try:
        # check columns
        actual_columns = set(df.columns)

        missing_columns = required_columns - actual_columns
        if missing_columns:
            raise ValueError(f"missing required columns: {missing_columns}")

        # Check dataframe has any data rows
        if df.empty:
            raise ValueError(f"df has no data rows")

        return True

    except Exception as e:
        raise SystemExit(f"validation failed: {e}")

def boundaries_to_arrow_old(df_in: pd.DataFrame, out_dir: str = None) -> None:
    out_dir = Path(out_dir) / "arrow" / 'arrow_boundaries'
    out_dir.mkdir(parents=True, exist_ok=True)


    shards = []
    total_polys = 0
    total_points = 0

    for boundaries in df_in:
        x_lists: List[List[float]] = []
        y_lists: List[List[float]] = []
        plane_ids: List[int] = []
        cell_ids: List[int] = []

        for _, row in boundaries.iterrows():
            # transform to Arrow-compatible arrays with list columns for coordinates. For example from
            # plane_id=15, cell_id=1001, coords=[[100,200], [105,200], [105,205]] into Arrow list columns where
            # one polygon becomes x_list=[100,105,105], y_list=[200,200,205], plane_id=15, cell_id=1001.

            pid = row["plane_id"]
            cid = row["cell_id"]
            coords = row["coords"]
            if not coords:
                raise ValueError
            xs = [float(x) for x, _ in coords]
            ys = [float(y) for _, y in coords]
            if not xs:
                raise ValueError
            x_lists.append(xs)
            y_lists.append(ys)
            plane_ids.append(pid)
            cell_ids.append(cid)

        # Write one feather per plane file
        plane_suffix = f"{plane_ids[0]:02d}" if plane_ids else "00"
        arrays = {
            "x_list": pa.array(x_lists, type=pa.list_(pa.float32())),
            "y_list": pa.array(y_lists, type=pa.list_(pa.float32())),
            "plane_id": pa.array(pd.Series(plane_ids, dtype="uint16")),
            "label": pa.array(pd.Series(cell_ids, dtype="int32")),
        }
        table = pa.table(arrays)
        shard_name = f"boundaries_plane_{plane_suffix}.feather"
        feather.write_feather(table, (out_dir / shard_name).as_posix(), compression="uncompressed")
        polys = len(x_lists)
        pts = sum(len(xs) for xs in x_lists)
        total_polys += polys
        total_points += pts
        shards.append({"url": shard_name, "rows": int(polys), "plane": int(plane_ids[0] if plane_ids else -1)})
        # print(f"Wrote {shard_name}: polys={polys}, points={pts}")

    # Manifest
    manifest = {
        "format": "arrow-feather",
        "total_rows": int(total_polys),  # polygons count
        "total_points": int(total_points),
        "shards": shards,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    io_utils_logger.info(f"Done. Total polys: {total_polys}. Total points: {total_points}. Files: {len(shards)}. Output: {out_dir}")



def _boundaries_to_arrow(df_in: List[pd.DataFrame], out_dir: str = None) -> None:
    out_dir = Path(out_dir) / "arrow" / 'arrow_boundaries'
    out_dir.mkdir(parents=True, exist_ok=True)

    shards = []
    total_polys = 0
    total_points = 0

    for boundaries in df_in:
        if boundaries.empty:
            continue

        # Vectorized operations using NumPy arrays
        coords_array = boundaries["coords"].values
        plane_ids_array = boundaries["plane_id"].values
        cell_ids_array = boundaries["cell_id"].values

        # Pre-allocate lists with known size
        num_rows = len(boundaries)
        x_lists = []
        y_lists = []
        plane_ids = np.empty(num_rows, dtype=np.uint16)
        cell_ids = np.empty(num_rows, dtype=np.int32)

        # Vectorized coordinate extraction
        valid_idx = 0
        for i in range(num_rows):
            coords = coords_array[i]
            if not coords:  # Skip empty coordinates
                continue

            # Use NumPy for faster array operations
            coords_np = np.array(coords, dtype=np.float32)
            xs = coords_np[:, 0]
            ys = coords_np[:, 1]

            if len(xs) == 0:
                continue

            x_lists.append(xs.tolist())
            y_lists.append(ys.tolist())
            plane_ids[valid_idx] = plane_ids_array[i]
            cell_ids[valid_idx] = cell_ids_array[i]
            valid_idx += 1

        # Trim arrays to actual size
        plane_ids = plane_ids[:valid_idx]
        cell_ids = cell_ids[:valid_idx]

        if valid_idx == 0:
            continue

        # Create Arrow arrays directly without intermediate pandas Series
        plane_suffix = f"{plane_ids[0]:02d}"
        arrays = {
            "x_list": pa.array(x_lists, type=pa.list_(pa.float32())),
            "y_list": pa.array(y_lists, type=pa.list_(pa.float32())),
            "plane_id": pa.array(plane_ids),
            "label": pa.array(cell_ids),
        }

        table = pa.table(arrays)
        shard_name = f"boundaries_plane_{plane_suffix}.feather"

        # Write with optimal compression settings
        feather.write_feather(
            table,
            (out_dir / shard_name).as_posix(),
            compression="uncompressed"
        )

        # Calculate stats efficiently
        polys = len(x_lists)
        pts = sum(len(xs) for xs in x_lists)  # This is still the fastest way
        total_polys += polys
        total_points += pts

        shards.append({
            "url": shard_name,
            "rows": polys,
            "plane": int(plane_ids[0])
        })
        # print(f"Wrote {shard_name}: polys={polys}, points={pts}")

    # Manifest
    manifest = {
        "format": "arrow-feather",
        "total_rows": total_polys,  # No need for int() conversion
        "total_points": total_points,
        "shards": shards,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    io_utils_logger.info(f"Done. Total polys: {total_polys}. Total points: {total_points}. Files: {len(shards)}. Output: {out_dir}")



def boundaries_to_arrow(dfs_in: List[pd.DataFrame], out_dir: str, compression: str = "uncompressed"):
    """
    Converts a list of DataFrames of boundary data into one Arrow Feather file per plane.

    Args:
        dfs_in: A list of DataFrames. Each DataFrame must contain data for a single plane
                and have the columns ['plane_id', 'label', 'coords']. The 'coords' column
                should contain lists of [x, y] coordinates.
        out_dir: The root directory to save the output 'arrow_boundaries' folder to.
        compression: The compression to use for the Feather files.
    """
    out_dir = Path(out_dir) / "arrow" / 'arrow_boundaries'
    out_dir.mkdir(parents=True, exist_ok=True)

    comp = compression if compression != "none" else None

    shards = []
    total_polys = 0
    total_points = 0

    # Define the schema once, to be used for all files
    schema = pa.schema([
        pa.field('x_list', pa.list_(pa.float32())),
        pa.field('y_list', pa.list_(pa.float32())),
        pa.field('plane_id', pa.uint16()),
        pa.field('label', pa.int32())
    ])

    # Process each DataFrame in the input list
    for idx, df_plane in enumerate(dfs_in):
        # Check for required columns
        required_cols = ['plane_id', 'cell_id', 'coords']
        if not all(col in df_plane.columns for col in required_cols):
            io_utils_logger.info("Warning: A DataFrame is missing required columns. Skipping.")
            continue

        # Get the plane ID - use index as fallback for empty DataFrames
        if df_plane.empty:
            # Empty plane - use the list index as plane_id
            current_plane_id = idx
        else:
            # Non-empty plane - get plane_id from first row
            current_plane_id = int(df_plane['plane_id'].iloc[0])

            # Validate: plane_id from data should match list index
            if current_plane_id != idx:
                raise ValueError(
                    f"Plane ID mismatch: DataFrame at index {idx} has plane_id={current_plane_id}. "
                    f"Expected plane_id to match index. Check that dfs_in is ordered correctly by plane."
                )

        shard_name = f"boundaries_plane_{current_plane_id:02d}.feather"

        # Filter out rows with empty coordinate lists (skip if already empty)
        if not df_plane.empty:
            df_plane = df_plane.copy()
            df_plane = df_plane[df_plane["coords"].str.len() > 0]

        if df_plane.empty:
            # If plane has no valid polygons, write an empty Feather file
            empty_table = schema.empty_table()
            feather.write_feather(empty_table, (out_dir / shard_name).as_posix(), compression=comp)
            shards.append({"url": shard_name, "rows": 0, "plane": current_plane_id})
            # io_utils_logger.info(f"Wrote empty shard {shard_name} for plane {current_plane_id}")
            continue

        # Prepare data for Arrow, using the 'coords' column directly
        x_lists = df_plane["coords"].apply(lambda coords: [float(x) for x, _ in coords])
        y_lists = df_plane["coords"].apply(lambda coords: [float(y) for _, y in coords])
        labels = pd.to_numeric(df_plane["cell_id"], errors="coerce").fillna(-1).astype("int32")

        # Create the Arrow table
        arrays = {
            "x_list": pa.array(x_lists.tolist(), type=pa.list_(pa.float32())),
            "y_list": pa.array(y_lists.tolist(), type=pa.list_(pa.float32())),
            "plane_id": pa.array([current_plane_id] * len(df_plane), type=pa.uint16()),
            "label": pa.array(labels.tolist(), type=pa.int32()),
        }
        table = pa.table(arrays, schema=schema)

        # Write the Feather file
        feather.write_feather(table, (out_dir / shard_name).as_posix(), compression=comp)

        polys = len(df_plane)
        pts = sum(x_lists.str.len())
        total_polys += polys
        total_points += pts

        shards.append({"url": shard_name, "rows": int(polys), "plane": current_plane_id})
        # print(f"Wrote {shard_name}: polys={polys}, points={pts}")

    # Manifest
    manifest = {
        "format": "arrow-feather",
        "total_rows": int(total_polys),
        "total_points": int(total_points),
        "shards": sorted(shards, key=lambda s: s['plane']),  # Sort shards by plane number
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    # io_utils_logger.info(f"Done. Total polys: {total_polys}. Total points: {total_points}. Files: {len(shards)}. Output: {out_dir}")
    io_utils_logger.info(f"Saved at: {out_dir}")









