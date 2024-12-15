from typing import Tuple, Optional, Dict, Any, Union
import logging
import os
import pickle
import sys
import tempfile
import redis
from pathlib import Path
import shutil
from urllib.parse import urlparse
from urllib.request import urlopen

import numpy as np
import pandas as pd
import numpy_groupies as npg
from tqdm import tqdm

from pciSeq.src.diagnostics.utils import RedisDB
from pciSeq.src.preprocess.utils import get_img_shape
from pciSeq.src.viewer.utils import copy_viewer_code, make_config_js, make_classConfig_js

utils_logger = logging.getLogger(__name__)


def log_file(cfg):
    """
    Setup the logger file handler if it doesn't already exist.

    Args:
        cfg: Configuration dictionary containing 'output_path' key

    Notes:
        - Creates a FileHandler if one doesn't exist
        - Writes logs to 'pciSeq.log' in the output directory
        - Uses format: '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    """
    root_logger = logging.getLogger()
    if root_logger.handlers:
        # setup a FileHandler if it has not been setup already. Maybe I should be adding a FileHandler anyway,
        # regardless whether there is one already or not
        if not np.any([isinstance(d, logging.FileHandler) for d in root_logger.handlers]):
            # Add type check for output_path
            if not isinstance(cfg['output_path'], str):
                raise TypeError("cfg['output_path'] must be a string")
            logfile = os.path.join(get_out_dir(cfg['output_path']), 'pciSeq.log')
            fh = logging.FileHandler(logfile, mode='w')
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            fh.setFormatter(formatter)

            root_logger.addHandler(fh)
            utils_logger.info('Writing to %s' % logfile)


def check_redis_server() -> bool:
    """
    Check if Redis server is running and available.

    Returns:
        bool: True if Redis server is running, False otherwise

    Notes:
        - Attempts to create a RedisDB connection
        - Logs warning if connection fails
        - Required for diagnostics functionality

    Example:
        >>> if check_redis_server():
        ...     print("Redis server is running")
        ... else:
        ...     print("Redis server is not available")
    """
    utils_logger.info("Checking Redis server connection...")
    try:
        RedisDB()
        utils_logger.info("Redis server connection successful")
        return True
    except (redis.exceptions.ConnectionError, ConnectionRefusedError, OSError) as e:
        utils_logger.warning(
            "Redis server not available. Diagnostics will be disabled. "
            f"Error: {str(e)}"
        )
        return False


def write_data(
        cell_data: pd.DataFrame,
        gene_data: pd.DataFrame,
        cell_boundaries: pd.DataFrame,
        var_bayes: Any,
        cfg: Dict[str, Any]
) -> None:
    """
    Write analysis results to output files.

    Args:
        cell_data: DataFrame containing cell information
        gene_data: DataFrame containing gene information
        cell_boundaries: DataFrame containing cell boundary coordinates
        var_bayes: Variational Bayes results object
        cfg: Configuration dictionary with 'output_path' and 'InsideCellBonus'

    Notes:
        - Creates 'data' subdirectory in output path
        - Saves files in TSV format
        - Serializes var_bayes object for debugging

    Files created:
        - cellData.tsv: Cell information
        - geneData.tsv: Gene information
        - cellBoundaries.tsv: Cell boundary coordinates
        - debug/pciSeq.pickle: Serialized var_bayes object
    """
    try:
        # Create output directory
        out_dir = Path(get_out_dir(cfg['output_path'])) / 'data'
        out_dir.mkdir(parents=True, exist_ok=True)

        # Save cell data
        cell_path = out_dir / 'cellData.tsv'
        cell_data.to_csv(cell_path, sep='\t', index=False)
        utils_logger.info(f'Saved cell data to {cell_path}')

        # Save gene data
        gene_path = out_dir / 'geneData.tsv'
        gene_data.to_csv(gene_path, sep='\t', index=False)
        utils_logger.info(f'Saved gene data to {gene_path}')

        # Save boundaries based on configuration
        bounds_path = out_dir / 'cellBoundaries.tsv'
        if not cfg['InsideCellBonus']:
            # Use gaussian contours
            ellipsoid_boundaries = cell_data[['Cell_Num', 'gaussian_contour']]
            ellipsoid_boundaries = ellipsoid_boundaries.rename(
                columns={"Cell_Num": "cell_id", "gaussian_contour": "coords"}
            )
            ellipsoid_boundaries.to_csv(bounds_path, sep='\t', index=False)
        else:
            # Use segmentation borders
            cell_boundaries.to_csv(bounds_path, sep='\t', index=False)
        utils_logger.info(f'Saved cell boundaries to {bounds_path}')

        # Serialize var_bayes object
        serialise(var_bayes, out_dir / 'debug')

    except Exception as e:
        utils_logger.error(f"Failed to write data: {str(e)}")
        raise


def serialise(var_bayes: Any, debug_dir: Union[str, Path]) -> None:
    """
    Serialize variational Bayes object for debugging purposes.

    Args:
        var_bayes: Variational Bayes results object to serialize
        debug_dir: Directory path to save the pickle file

    Notes:
        - Creates debug directory if it doesn't exist
        - Saves object as 'pciSeq.pickle'
    """
    debug_dir = Path(debug_dir)
    debug_dir.mkdir(parents=True, exist_ok=True)

    pickle_path = debug_dir / 'pciSeq.pickle'
    try:
        with open(pickle_path, 'wb') as outf:
            pickle.dump(var_bayes, outf)
            utils_logger.info(f'Serialized object saved to {pickle_path}')
    except Exception as e:
        utils_logger.error(f"Failed to serialize object: {str(e)}")
        raise


def pre_launch(
        cell_data: pd.DataFrame,
        coo: Any,
        scRNA_seq: Optional[Any],
        cfg: Dict[str, Any]
) -> Path:
    """
    Perform pre-flight checks and setup before launching the viewer.

    Args:
        cell_data: DataFrame containing cell information
        coo: Coordinate data object
        scRNA_seq: Single-cell RNA sequencing data (optional)
        cfg: Configuration dictionary

    Returns:
        Path: Destination folder containing viewer code

    Notes:
        - Gets image dimensions from coordinates
        - Copies viewer code to destination
        - Creates configuration JavaScript files
    """
    try:
        # Get image dimensions
        height, width = get_img_shape(coo)

        # Setup destination directory
        dst = Path(get_out_dir(cfg['output_path']))

        # Copy viewer code and create configs
        copy_viewer_code(cfg, dst)
        make_config_js(dst, width, height)

        # Handle class configuration
        if scRNA_seq is None:
            # Extract unique class names
            label_list = [d[:] for d in cell_data.ClassName.values]
            labels = sorted(set(item for sublist in label_list for item in sublist))
            make_classConfig_js(labels, dst)

        return dst

    except Exception as e:
        utils_logger.error(f"Pre-launch setup failed: {str(e)}")
        raise


def negative_binomial_loglikelihood(x: np.ndarray, r: float, p: np.ndarray) -> np.ndarray:
    """
    Calculate the Negative Binomial log-likelihood for given parameters.

    The Negative Binomial distribution models the number of successes (x) before
    r failures occur, with probability of success p. The PMF is:
    P(X=k) = C(k+r-1,k) * p^k * (1-p)^r

    Args:
        x: Array of observed counts with shape
        r: Number of failures until stopping (dispersion parameter)
        p: Probability of success, array with shape

    Returns:
        Array of log-likelihood contributions with shape
    """
    try:
        x = x[:, :, None]  # Add dimension for broadcasting

        # negative binomial log-likelihood.
        # Focusing only on the terms that involve p and r (without the binomial coefficient):
        log_likelihood = x * np.log(p) + r * np.log(1 - p)

        return log_likelihood

    except Exception as e:
        utils_logger.error(f"Error calculating negative binomial log-likelihood: {str(e)}")
        raise ValueError("Failed to compute log-likelihood. Check input dimensions and values.")


def has_converged(
        spots: Any,
        p0: Optional[np.ndarray],
        tol: float
) -> Tuple[bool, float]:
    """
    Check if the probability assignments have converged.

    Args:
        spots: Spot data object containing parent_cell_prob
        p0: Previous probability matrix (None for first iteration)
        tol: Convergence tolerance threshold

    Returns:
        Tuple containing:
            - bool: True if converged, False otherwise
            - float: Maximum absolute difference between iterations
    """
    p1 = spots.parent_cell_prob
    if p0 is None:
        p0 = np.zeros_like(p1)

    try:
        delta = np.max(np.abs(p1 - p0))
        converged = (delta < tol)
        return converged, delta
    except Exception as e:
        utils_logger.error(f"Convergence check failed: {str(e)}")
        raise


def download_url_to_file(
        url: str,
        dst: Union[str, Path],
        progress: bool = True
) -> None:
    """
    Download a file from URL to local path with progress tracking.

    Args:
        url: URL of the file to download
        dst: Destination path for downloaded file
        progress: Whether to display progress bar (default: True)

    Notes:
        - Uses temporary file for safe downloading
        - Shows progress with tqdm if enabled
        - Cleans up temporary files on failure
    """
    BUFFER_SIZE = 8192

    dst = Path(dst).expanduser()
    dst_dir = dst.parent

    try:
        # Get file size if available
        with urlopen(url) as u:
            meta = u.info()
            file_size = int(meta.get("Content-Length", 0))

        # Download with temporary file
        with tempfile.NamedTemporaryFile(delete=False, dir=dst_dir) as temp_file:
            with tqdm(
                    total=file_size,
                    disable=not progress,
                    unit='B',
                    unit_scale=True,
                    unit_divisor=1024
            ) as pbar:
                with urlopen(url) as response:
                    while True:
                        buffer = response.read(BUFFER_SIZE)
                        if not buffer:
                            break
                        temp_file.write(buffer)
                        pbar.update(len(buffer))

            # Move temporary file to destination
            shutil.move(temp_file.name, dst)
            utils_logger.info(f"Downloaded {url} to {dst}")

    except Exception as e:
        utils_logger.error(f"Download failed: {str(e)}")
        if 'temp_file' in locals():
            try:
                os.remove(temp_file.name)
            except OSError:
                pass
        raise


def load_from_url(url: str) -> str:
    """
    Download and load a file from URL if not already present.

    Args:
        url: URL of the file to download

    Returns:
        str: Path to the downloaded file
    """
    try:
        parts = urlparse(url)
        filename = Path(parts.path).name

        if not Path(filename).exists():
            sys.stderr.write(f'Downloading: "{url}" to {filename}\n')
            download_url_to_file(url, filename)

        return filename

    except Exception as e:
        utils_logger.error(f"Failed to load from URL: {str(e)}")
        raise


def get_out_dir(path: Optional[str] = None, sub_folder: str = '') -> str:
    """
    Get or create output directory path.

    Args:
        path: List containing base path, or None for default temp directory
        sub_folder: Optional subdirectory name

    Returns:
        str: Path to output directory

    Notes:
        - Uses system temp directory if path is None or ['default']
        - Creates directories if they don't exist
    """
    if path is None or path == 'default':
        out_dir = Path(tempfile.gettempdir()) / 'pciSeq'
    else:
        out_dir = Path(path) / sub_folder / 'pciSeq'

    out_dir.mkdir(parents=True, exist_ok=True)
    return str(out_dir)


def gaussian_contour(
        mu: np.ndarray,
        cov: np.ndarray,
        sdwidth: float = 3.0
) -> np.ndarray:
    """
    Generate ellipsoid contour points from covariance matrix and mean.

    Args:
        mu: Mean vector [2,] specifying the center of the ellipse
        cov: Covariance matrix [2,2] specifying the shape and orientation
        sdwidth: Number of standard deviations for contour size (default: 3.0)

    Returns:
        Array of contour points shape [n_points, 2]

    Raises:
        ValueError: If input dimensions are incorrect

    Examples:
        Basic usage with circular covariance:
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> mu = np.array([0, 0])
        >>> cov = np.array([[1, 0], [0, 1]])
        >>> points = gaussian_contour(mu, cov)
        >>> plt.plot(points[:, 0], points[:, 1])

        Elliptical shape with correlation:
        >>> mu = np.array([1, 2])
        >>> cov = np.array([[2, 0.5],
        ...                 [0.5, 1]])
        >>> points = gaussian_contour(mu, cov)
        >>> plt.plot(points[:, 0], points[:, 1])

        Multiple ellipses with different parameters:
        >>> fig, ax = plt.subplots()
        >>> # Standard circle
        >>> p1 = gaussian_contour([0, 0], [[1, 0], [0, 1]])
        >>> # Elongated ellipse
        >>> p2 = gaussian_contour([2, 2], [[3, 0.5], [0.5, 1]])
        >>> # Rotated ellipse
        >>> p3 = gaussian_contour([-1, 1], [[1, -0.8], [-0.8, 1]])
        >>> ax.plot(p1[:, 0], p1[:, 1], 'b-', label='Circle')
        >>> ax.plot(p2[:, 0], p2[:, 1], 'r-', label='Elongated')
        >>> ax.plot(p3[:, 0], p3[:, 1], 'g-', label='Rotated')
        >>> ax.legend()
        >>> ax.axis('equal')
        >>> plt.show()

    Notes:
        - The covariance matrix must be 2x2 symmetric and positive semi-definite
        - sdwidth controls the size of the ellipse (larger values = larger ellipse)
        - The contour represents a constant probability density contour of a
          2D Gaussian distribution
    """
    if mu.shape != (2,) or cov.shape != (2, 2):
        raise ValueError("Invalid input dimensions. Expected mu.shape=(2,) and cov.shape=(2,2)")

    try:
        contour_points = 40

        # Generate circle points
        tt = np.linspace(0, 2 * np.pi, contour_points)
        ap = np.stack([np.cos(tt), np.sin(tt)])

        # Transform to ellipse
        eigvals, eigvecs = np.linalg.eigh(cov)
        eigvals = sdwidth * np.sqrt(np.maximum(eigvals, 0))  # Ensure non-negative
        vd = eigvecs @ (eigvals[:, None] * ap)
        out = vd + mu.reshape(2, -1)

        return np.array(list(zip(*out)), dtype=np.float32)
    except Exception as e:
        utils_logger.error(f"Error generating gaussian contour: {str(e)}")
        raise


def read_image_objects(img_obj, cfg):
    meanCellRadius = np.mean(np.sqrt(img_obj.area / np.pi)) * 0.5
    relCellRadius = np.sqrt(img_obj.area / np.pi) / meanCellRadius

    # append 1 for the misreads
    relCellRadius = np.append(1, relCellRadius)

    InsideCellBonus = cfg['InsideCellBonus']
    if not InsideCellBonus:
        # This is more for clarity. The operation below will work fine even if InsideCellBonus is False
        InsideCellBonus = 0

    # if InsideCellBonus == 0 then CellAreaFactor will be equal to 1.0
    numer = np.exp(-relCellRadius ** 2 / 2) * (1 - np.exp(InsideCellBonus)) + np.exp(InsideCellBonus)
    denom = np.exp(-0.5) * (1 - np.exp(InsideCellBonus)) + np.exp(InsideCellBonus)
    CellAreaFactor = numer / denom

    out = {}
    out['area_factor'] = CellAreaFactor.astype(np.float32)
    out['rel_radius'] = relCellRadius.astype(np.float32)
    out['area'] = np.append(np.nan, img_obj.area.astype(np.uint32))
    out['x0'] = np.append(-sys.maxsize, img_obj.x0.values).astype(np.float32)
    out['y0'] = np.append(-sys.maxsize, img_obj.y0.values).astype(np.float32)
    out['cell_label'] = np.append(0, img_obj.label.values).astype(np.uint32)
    if 'old_label' in img_obj.columns:
        out['cell_label_old'] = np.append(0, img_obj.old_label.values).astype(np.uint32)
    # First cell is a dummy cell, a super neighbour (ie always a neighbour to any given cell)
    # and will be used to get all the misreads. It was given the label=0 and some very small
    # negative coords

    return out, meanCellRadius.astype(np.float32)

def keep_labels_unique(scdata):
    """
    In the single cell data you might find cases where two or more rows have the same gene label
    In these cases keep the row with the highest total gene count
    """

    # 1. get the row total and assign it to a new column
    scdata = scdata.assign(total=scdata.sum(axis=1))

    # 2. rank by gene label and total gene count and keep the one with the highest total
    scdata = scdata.sort_values(['gene_name', 'total'], ascending=[True, False]).groupby('gene_name').head(1)

    # 3. Drop the total column and return
    return scdata.drop(['total'], axis=1)


# @dask.delayed
def scaled_exp(
        cell_area_factor: np.ndarray,
        sc_mean_expressions: np.ndarray,
        inefficiency: np.ndarray
) -> np.ndarray:
    """
    Calculate scaled expression values using cell areas and efficiency factors.

    Args:
        cell_area_factor: Array of cell area scaling factors
        sc_mean_expressions: Array of mean expression values
        inefficiency: Array of inefficiency factors

    Returns:
        Array of scaled expression values

    Optimization:
        - Use np.einsum for efficient matrix multiplication
        - Pre-allocate output array for large computations
        - Consider using numba for very large arrays
        - Add input validation for array shapes
    """
    try:
        if np.all(cell_area_factor == 1):
            subscripts = 'gk,g->gk'
            operands = [sc_mean_expressions, inefficiency]
        else:
            subscripts = 'c,gk,g->cgk'
            operands = [cell_area_factor, sc_mean_expressions, inefficiency]

        return np.einsum(subscripts, *operands)

    except Exception as e:
        utils_logger.error(f"Error in scaled expression calculation: {str(e)}")
        raise


def index_genes(gene_names: np.ndarray):
    """Creates a mapping between gene names and numerical indices.

    This function serves two purposes:
    1. Returns unique gene names (used by Genes class to establish gene panel)
    2. Maps each spot to a gene_id
    3. the gene_id is the position of the gene reading name in the ordered gene names array

    Args:
        gene_names: array of gene names

    Returns:
        - unique_genes: sorted array of unique gene names
        - gene_ids: numerical index for each original gene name (maps each spot to its gene type)
        - counts: frequency of each unique gene
    """
    return np.unique(gene_names, return_inverse=True, return_counts=True)


def empirical_mean(spots, cells):
    # spots = self.spots

    # get the total gene counts per cell
    N_c = cells.total_counts

    xy_spots = spots.xy_coords
    prob = spots.parent_cell_prob
    n = cells.config['nNeighbors'] + 1

    # multiply the x coord of the spots by the cell prob
    a = np.tile(xy_spots[:, 0], (n, 1)).T * prob

    # multiply the y coord of the spots by the cell prob
    b = np.tile(xy_spots[:, 1], (n, 1)).T * prob

    # aggregated x and y coordinate
    idx = spots.parent_cell_id
    x_agg = npg.aggregate(idx.ravel(), a.ravel(), size=len(N_c))
    y_agg = npg.aggregate(idx.ravel(), b.ravel(), size=len(N_c))

    # get the estimated cell centers
    x_bar = np.nan * np.ones(N_c.shape)
    y_bar = np.nan * np.ones(N_c.shape)

    x_bar[N_c > 0] = x_agg[N_c > 0] / N_c[N_c > 0]
    y_bar[N_c > 0] = y_agg[N_c > 0] / N_c[N_c > 0]

    # cells with N_c = 0 will end up with x_bar = y_bar = np.nan
    xy_bar_fitted = np.array(list(zip(x_bar.T, y_bar.T)))

    # if you have a value for the estimated centroid use that, otherwise
    # use the initial (starting values) centroids
    ini_centr = cells.ini_centroids()
    xy_bar = np.array(tuple(zip(*[ini_centr['x'], ini_centr['y']])))

    # # sanity check. NaNs or Infs should appear together
    # assert np.all(np.isfinite(x_bar) == np.isfinite(y_bar))
    # use the fitted centroids where possible otherwise use the initial ones
    xy_bar[np.isfinite(x_bar)] = xy_bar_fitted[np.isfinite(x_bar)]
    # self.cells.centroid = pd.DataFrame(xy_bar, columns=['x', 'y'])
    return pd.DataFrame(xy_bar, columns=['x', 'y'])