"""File and data I/O operations."""
import os
import numpy as np
from pathlib import Path
import pickle
import tempfile
import shutil
import sqlite3
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
logger = logging.getLogger(__name__)


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
        logger.info('Writing to %s' % logfile)


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
        logger.info('Downloading: "%s" to %s', url, filename)
        download_url_to_file(url, filename)
    return filename


def collect_metadata() -> Dict:
    """Collect metadata about the environment and analysis run."""
    import platform
    import subprocess
    import sys
    from datetime import datetime

    # Git commit of the pciSeq code
    git_commit = None
    try:
        pciSeq_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)
        ))))
        git_commit = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'],
            cwd=pciSeq_dir,
            stderr=subprocess.DEVNULL,
            text=True
        ).strip()
    except Exception:
        pass

    # Key package versions
    pkg_versions = {}
    for pkg in ['numpy', 'scipy', 'pandas', 'pciSeq']:
        try:
            mod = __import__(pkg)
            pkg_versions[pkg] = getattr(mod, '__version__', 'unknown')
        except ImportError:
            pass

    metadata = {
        'date': datetime.now().isoformat(),
        'git_commit': git_commit,
        'hostname': platform.node(),
        'os': f'{platform.system()} {platform.release()}',
        'python_version': sys.version.split()[0],
        'package_versions': pkg_versions,
    }
    return metadata


def serialise(varBayes: Any, debug_dir: str) -> None:
    """Pickle variable Bayes object to debug directory.

    Args:
        varBayes: Object to serialize
        debug_dir: Directory to save pickle file
    """
    varBayes._metadata = collect_metadata()

    if not os.path.exists(debug_dir):
        os.makedirs(debug_dir)
    pickle_dst = os.path.join(debug_dir, 'pciSeq.pickle')
    with open(pickle_dst, 'wb') as outf:
        pickle.dump(varBayes, outf)

    pickle_mb = os.path.getsize(pickle_dst) / (1024 * 1024)
    logger.info('Saved at %s (%.1f MB)', pickle_dst, pickle_mb)

    # Export diagnostics database to diagnostics folder (sibling of arrow folder)
    # This allows the viewer to auto-discover it alongside arrow data
    data_dir = os.path.dirname(debug_dir)  # Go up from debug to data folder
    export_diagnostics(varBayes, data_dir)


def export_diagnostics(varBayes: Any, output_dir: str) -> None:
    """Export diagnostics data (check_cell and check_spot) to a single SQLite database.

    Writes to {output_dir}/diagnostics/diagnostics.db

    Tables:
      - metadata: key-value pairs (including JSON arrays)
      - cells: per-cell diagnostic data
      - spots: per-spot diagnostic data

    Args:
        varBayes: Fitted VarBayes object
        output_dir: Base data directory
    """
    diagnostics_dir = os.path.join(output_dir, 'diagnostics')
    os.makedirs(diagnostics_dir, exist_ok=True)

    db_path = os.path.join(diagnostics_dir, 'diagnostics.db')
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # --- Create Tables ---
    cursor.execute('CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT)')

    cursor.execute('''
        CREATE TABLE cells (
            cell_id INTEGER PRIMARY KEY,
            scaled_means BLOB,
            theta_bar BLOB,
            gene_count BLOB,
            class_prob BLOB
        )
    ''')

    cursor.execute('''
        CREATE TABLE spots (
            spot_id INTEGER PRIMARY KEY,
            gene_idx INTEGER,
            x INTEGER,
            y INTEGER,
            z INTEGER,
            neighbor_cell_ids TEXT,
            mvn_loglik BLOB,
            attention BLOB,
            expr_fluct BLOB,
            cell_inefficiency BLOB
        )
    ''')

    # --- Gather Data ---
    cells = varBayes.cells
    genes = varBayes.genes
    spots = varBayes.spots

    # Label Map
    label_map = {}
    if varBayes.config.get('label_map'):
        label_map = {str(k): int(v) for k, v in varBayes.config['label_map'].items()}

    # Gene Panel
    gene_panel = genes.gene_panel.tolist()

    # Misread Density
    misread_series = genes.misread_density
    if hasattr(misread_series, 'to_dict'):
        misread_dict = {str(k): float(v) for k, v in misread_series.to_dict().items()}
    else:
        logger.error("Diagnostics export skipped: 'misread_density' is missing or invalid.")
        return

    # --- Populate Metadata ---
    # Compute scaled_means for metadata nC (and for cells table)
    # logger.info('Computing scaled_exp for diagnostics export...')
    scaled_means = varBayes.scaled_exp.compute()
    nC, nG, nK = scaled_means.shape

    nS = spots.nS
    # Check if we can get nN (needs neighbor_ids)
    neighbor_ids = spots.parent_cell_id
    nN = 0
    if neighbor_ids is not None:
        nN = neighbor_ids.shape[1]

    meta_items = [
        # Cell-related
        ('nC', str(nC)),
        ('nG', str(nG)),
        ('nK', str(nK)),
        ('rSpot', str(float(varBayes.config['rSpot']))),
        ('SpotReg', str(float(varBayes.config['SpotReg']))),
        ('class_names', json.dumps(cells.class_names.tolist())),
        ('eta_bar', json.dumps(genes.eta_bar.astype(np.float32).tolist())),
        ('mean_gene_reads_per_class', json.dumps(cells.mean_gene_reads_per_class().astype(np.float32).tolist())),

        # Spot-related
        ('nS', str(int(nS))),
        ('nN', str(int(nN))),
        ('misread_density', json.dumps(misread_dict)),

        # Shared
        ('gene_panel', json.dumps(gene_panel)),
        ('label_map', json.dumps(label_map)),
    ]
    cursor.executemany('INSERT INTO metadata VALUES (?, ?)', meta_items)
    # logger.info('Inserted %d metadata entries', len(meta_items))

    # --- Populate Cells Table ---
    scaled_means_f32 = scaled_means.astype(np.float32)
    theta_bar_f32 = cells.theta_bar.astype(np.float32)
    gene_count_f32 = cells.geneCount.astype(np.float32)
    class_prob_f32 = cells.classProb.astype(np.float32)

    batch_size = 10000
    for batch_start in range(0, nC, batch_size):
        batch_end = min(batch_start + batch_size, nC)
        batch_data = []
        for c in range(batch_start, batch_end):
            batch_data.append((
                c,
                scaled_means_f32[c].tobytes(),
                theta_bar_f32[c].tobytes(),
                gene_count_f32[c].tobytes(),
                class_prob_f32[c].tobytes(),
            ))
        cursor.executemany('INSERT INTO cells VALUES (?, ?, ?, ?, ?)', batch_data)
        # if (batch_end % 10000 == 0) or (batch_end == nC):
        #     logger.info('Inserted %d/%d cells', batch_end, nC)

    # --- Populate Spots Table ---
    if spots.mvn_loglik_arr is None or spots.attention is None or spots.expr_fluctuations is None or spots.cell_inefficiency is None or neighbor_ids is None:
        logger.warning('check_spot data missing; spots table will be empty.')
    else:
        mvn_f32 = spots.mvn_loglik_arr.astype(np.float32)
        attn_f32 = spots.attention.astype(np.float32)
        expr_f32 = spots.expr_fluctuations.astype(np.float32)
        cineff_f32 = spots.cell_inefficiency.astype(np.float32)
        gene_idx = spots.gene_id.astype(np.int32)
        xs = spots.data['x'].astype(np.int32).to_numpy()
        ys = spots.data['y'].astype(np.int32).to_numpy()
        zs = spots.data['z'].astype(np.int32).to_numpy()

        batch_size = 10000
        for start in range(0, nS, batch_size):
            end = min(start + batch_size, nS)
            batch = []
            for i in range(start, end):
                neigh_json = json.dumps(list(map(int, neighbor_ids[i].tolist())))
                batch.append((
                    int(spots.data.index[i]),
                    int(gene_idx[i]),
                    int(xs[i]), int(ys[i]), int(zs[i]),
                    neigh_json,
                    mvn_f32[i].tobytes(),
                    attn_f32[i].tobytes(),
                    expr_f32[i].tobytes(),
                    cineff_f32[i].tobytes(),
                ))
            cursor.executemany('''
                INSERT INTO spots (spot_id, gene_idx, x, y, z, neighbor_cell_ids, mvn_loglik, attention, expr_fluct, cell_inefficiency)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', batch)
            # if (end % 50000 == 0) or (end == nS):
            #     logger.info('Inserted %d/%d spots', end, nS)

    conn.commit()
    conn.close()

    db_size_mb = os.path.getsize(db_path) / (1024 * 1024)
    logger.info('Saved at: %s (%.1f MB)', db_path, db_size_mb)


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
    logger.info('Saved at %s', fname)


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
    cellData_path = os.path.join(out_dir, "cellData.tsv")
    cellData.to_csv(cellData_path, sep='\t', index=False)
    cellData_mb = os.path.getsize(cellData_path) / (1024 * 1024)
    logger.info('Saved at: %s (%.1f MB)', cellData_path, cellData_mb)

    # Save gene data
    geneData_path = os.path.join(out_dir, "geneData.tsv")
    geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
    geneData_mb = os.path.getsize(geneData_path) / (1024 * 1024)
    logger.info('Saved at: %s (%.1f MB)', geneData_path, geneData_mb)

    # Save boundaries
    cellBoundaries_path = os.path.join(out_dir, "cellBoundaries.tsv")
    cellBoundaries.to_csv(cellBoundaries_path, sep='\t', index=False)
    cellBoundaries_mb = os.path.getsize(cellBoundaries_path) / (1024 * 1024)
    logger.info('Saved at %s: (%.1f MB)', cellBoundaries_path, cellBoundaries_mb)


def write_arrow(geneData:pd.DataFrame, cellData:pd.DataFrame, cellBoundaries:pd.DataFrame, out_dir: str = None) -> None:
    geneData_to_arrow(geneData, out_dir)
    cellData_to_arrow(cellData, out_dir)
    # logger.info('boundaries_to_arrow_old - Starting')
    # boundaries_to_arrow_old(cellBoundaries, out_dir)
    # logger.info('boundaries_to_arrow_old - Ending')

    # logger.info('boundaries_to_arrow - Starting')
    boundaries_to_arrow(cellBoundaries, out_dir)
    # logger.info('boundaries_to_arrow - Ending')

    # logger.info('Saved at %s', os.path.join(out_dir, 'cellBoundaries.tsv'))


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

    # logger.info(f"Saved {total_rows} rows in {len(shards)} shards at {out_dir}")
    logger.info(f"Saved at {out_dir}")



def cellData_to_arrow(df_in: pd.DataFrame, out_dir: str = None) -> None:
    """
    Convert cell data DataFrame to Arrow Feather shards for JavaScript viewer.

    Output schema (8 columns):
      - cell_id: int32 (cell identifier)
      - X, Y, Z: float32 (cell centroid coordinates)
      - class_name: list<string> (predicted cell types, ordered by probability)
      - prob: list<float32> (classification probabilities)
      - gene_names: list<string> (detected genes, ordered by expression)
      - gene_counts: list<float32> (gene expression counts)

    Args:
        df_in: Cell data DataFrame from cells_summary()
        out_dir: Base output directory (arrow/arrow_cells/ will be appended)

    Raises:
        ValueError: If required source columns are missing
    """
    out_dir = Path(out_dir) / "arrow" / "arrow_cells"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Validate required source columns - fail fast if missing
    required_source_cols = ["Cell_Num", "X", "Y", "Z", "ClassName", "Prob", "Genenames", "CellGeneCount"]
    missing = [col for col in required_source_cols if col not in df_in.columns]
    if missing:
        raise ValueError(
            f"cellData_to_arrow: Missing required source columns: {missing}. "
            f"Cannot generate Arrow files for viewer. Check cells_summary() output."
        )

    # Output schema: source column -> (output name, arrow type)
    schema_map = [
        ("Cell_Num", "cell_id", pa.int32()),
        ("X", "X", pa.float32()),
        ("Y", "Y", pa.float32()),
        ("Z", "Z", pa.float32()),
        ("ClassName", "class_name", pa.list_(pa.string())),
        ("Prob", "prob", pa.list_(pa.float32())),
        ("Genenames", "gene_names", pa.list_(pa.string())),
        ("CellGeneCount", "gene_counts", pa.list_(pa.float32())),
    ]

    shards = []
    total_rows = 0
    shard_index = 0
    chunk_size = 100_000

    for start in range(0, len(df_in), chunk_size):
        df = df_in.iloc[start:start + chunk_size]
        n = len(df)

        # Build arrays in fixed column order
        arrays = {}
        for src_col, out_col, arrow_type in schema_map:
            if isinstance(arrow_type, pa.ListType):
                # List columns: convert DataFrame lists to Arrow list arrays
                arrays[out_col] = pa.array(df[src_col].tolist(), type=arrow_type)
            else:
                # Scalar columns: explicit type casting
                arrays[out_col] = pa.array(df[src_col], type=arrow_type)

        # Create table with fixed column order
        col_names = [out_col for _, out_col, _ in schema_map]
        table = pa.table({col: arrays[col] for col in col_names})

        # Write shard
        shard_name = f"cells_shard_{shard_index:03d}.feather"
        feather.write_feather(table, (out_dir / shard_name).as_posix(), compression="uncompressed")

        shards.append({"url": shard_name, "rows": int(n)})
        total_rows += n
        shard_index += 1

    # Write manifest
    manifest = {"format": "arrow-feather", "total_rows": int(total_rows), "shards": shards}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    # logger.info(f"Saved {total_rows} cell records in {len(shards)} shards at {out_dir}")
    logger.info(f"Saved at {out_dir}")


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
    logger.info(f"Done. Total polys: {total_polys}. Total points: {total_points}. Files: {len(shards)}. Output: {out_dir}")



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
    logger.info(f"Done. Total polys: {total_polys}. Total points: {total_points}. Files: {len(shards)}. Output: {out_dir}")



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
            logger.info("Warning: A DataFrame is missing required columns. Skipping.")
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
            # logger.info(f"Wrote empty shard {shard_name} for plane {current_plane_id}")
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
    # logger.info(f"Done. Total polys: {total_polys}. Total points: {total_points}. Files: {len(shards)}. Output: {out_dir}")
    logger.info(f"Saved at: {out_dir}")







