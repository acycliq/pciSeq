"""File and data I/O operations."""
import os
from pathlib import Path
import pickle
import tempfile
import shutil
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
               cellBoundaries: pd.DataFrame, varBayes: Any, cfg: Dict) -> None:
    """Write all data files to output directory.

    Args:
        cellData: Cell data DataFrame
        geneData: Gene data DataFrame
        cellBoundaries: Cell boundaries DataFrame
        varBayes: Variable Bayes object
        cfg: Configuration dictionary
    """
    dst = get_out_dir(cfg['output_path'])
    out_dir = os.path.join(dst, 'data')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Save cell data
    cellData.to_csv(os.path.join(out_dir, 'cellData.tsv'), sep='\t', index=False)
    io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'cellData.tsv'))

    # Save gene data
    geneData.to_csv(os.path.join(out_dir, 'geneData.tsv'), sep='\t', index=False)
    io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'geneData.tsv'))

    # Save boundaries based on InsideCellBonus setting
    ellipsoidBoundaries = cellData[['Cell_Num', 'gaussian_contour']]
    ellipsoidBoundaries = ellipsoidBoundaries.rename(
        columns={"Cell_Num": "cell_id", "gaussian_contour": "coords"})
    ellipsoidBoundaries.to_csv(os.path.join(out_dir, 'ellipsoidBoundaries.tsv'), sep='\t', index=False)

    cellBoundaries.to_csv(os.path.join(out_dir, 'cellBoundaries.tsv'), sep='\t', index=False)
    io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'cellBoundaries.tsv'))


    # flatfiles relevant to the cells that appear in the cellBoundaries tsv
    cellData = cellData.merge(cellBoundaries, left_on='Cell_Num', right_on='cell_id')
    cell_num_set = set(cellData['Cell_Num'])
    filtered_geneData = geneData[geneData['neighbour_array'].apply(lambda x: any(neighbour in cell_num_set for neighbour in x))]

    filtered_geneData.to_csv(os.path.join(out_dir, 'geneData_filtered.tsv'), sep='\t', index=False)
    io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'geneData_filtered.tsv'))

    cellData.iloc[:, -2:].to_csv(os.path.join(out_dir, 'cellBoundaries_filtered.tsv'), sep='\t', index=False)
    io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'cellBoundaries_filtered.tsv'))

    cellData.iloc[:, :-2].to_csv(os.path.join(out_dir, 'cellData_filtered.tsv'), sep='\t', index=False)
    io_utils_logger.info('Saved at %s', os.path.join(out_dir, 'cellData_filtered.tsv'))


    # Save debug info
    serialise(varBayes, os.path.join(out_dir, 'debug'))


