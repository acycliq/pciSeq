"""Statistical calculation utilities."""
import numpy as np
import pandas as pd
import numpy_groupies as npg
from typing import Tuple, Optional, Any, Union
import logging

# Configure logging
ops_utils_logger = logging.getLogger(__name__)


def negative_binomial_loglikelihood(x: np.ndarray, r: float, p: np.ndarray) -> np.ndarray:
    """Calculate the Negative Binomial log-likelihood for given parameters.

    The Negative Binomial distribution models the number of successes (x) before
    r failures occur, with probability of success p. The PMF is:
    P(X=k) = C(k+r-1,k) * p^k * (1-p)^r

    Args:
        x: Array of observed counts
        r: Number of failures until stopping (dispersion parameter)
        p: Probability of success array

    Returns:
        Array of log-likelihood contributions

    Raises:
        ValueError: If computation fails due to dimension mismatch or invalid values
    """
    try:
        x = x[:, :, None]  # Add dimension for broadcasting

        # negative binomial log-likelihood.
        # Focusing only on the terms that involve p and r (without the binomial coefficient):
        log_likelihood = x * np.log(p) + r * np.log(1 - p)

        return log_likelihood

    except Exception as e:
        ops_utils_logger.error(f"Error calculating negative binomial log-likelihood: {str(e)}")
        raise ValueError("Failed to compute log-likelihood. Check input dimensions and values.")


def softmax(X: np.ndarray, theta: float = 1.0, axis: Optional[int] = None) -> np.ndarray:
    """Compute the softmax of each element along an axis of X.

    Args:
        X: Input array (should be floats)
        theta: Multiplier prior to exponentiation (default: 1.0)
        axis: Axis to compute values along (default: first non-singleton axis)

    Returns:
        Array same size as X, normalized along the specified axis

    Notes:
        From https://nolanbconaway.github.io/blog/2017/softmax-numpy
    """
    # Make X at least 2d
    y = np.atleast_2d(X)

    # Find axis if not specified
    if axis is None:
        axis = next(j[0] for j in enumerate(y.shape) if j[1] > 1)

    # Multiply y against the theta parameter
    y = y * float(theta)

    # Subtract the max for numerical stability
    y = y - np.expand_dims(np.max(y, axis=axis), axis)

    # Exponentiate y
    y = np.exp(y)

    # Take the sum along the specified axis
    ax_sum = np.expand_dims(np.sum(y, axis=axis), axis)

    # Finally: divide elementwise
    p = y / ax_sum

    # Flatten if X was 1D
    if len(X.shape) == 1:
        p = p.flatten()

    return p


def has_converged(
        spots: Any,
        p0: Optional[np.ndarray],
        tol: float
) -> Tuple[bool, float]:
    """Check if probability assignments have converged.

    Args:
        spots: Spot data object containing parent_cell_prob
        p0: Previous probability matrix (None for first iteration)
        tol: Convergence tolerance threshold

    Returns:
        Tuple containing:
            - bool: True if converged, False otherwise
            - float: Maximum absolute difference between iterations

    Raises:
        Exception: If convergence check fails
    """
    p1 = spots.parent_cell_prob
    if p0 is None:
        p0 = np.zeros_like(p1)

    try:
        delta = np.max(np.abs(p1 - p0))
        converged = (delta < tol)
        return converged, delta
    except Exception as e:
        ops_utils_logger.error(f"Convergence check failed: {str(e)}")
        raise


def scaled_exp(cell_area_factor: np.ndarray,
               sc_mean_expressions: np.ndarray,
               inefficiency: np.ndarray) -> np.ndarray:
    """Calculate scaled expression values.

    Args:
        cell_area_factor: Cell area scaling factors
        sc_mean_expressions: Single cell mean expression values
        inefficiency: Inefficiency factors

    Returns:
        Scaled expression array
    """
    if np.all(cell_area_factor == 1):
        subscripts = 'gk, g -> gk'
        operands = [sc_mean_expressions, inefficiency]
    else:
        subscripts = 'c, gk, g -> cgk'
        operands = [cell_area_factor, sc_mean_expressions, inefficiency]

    return np.einsum(subscripts, *operands)


def empirical_mean(spots, cells):

    # get the total gene counts per cell
    N_c = cells.total_counts

    xyz_spots = spots.xyz_coords
    prob = spots.parent_cell_prob
    n = cells.config['nNeighbors'] + 1

    # multiply the x coord of the spots by the cell prob
    a = np.tile(xyz_spots[:, 0], (n, 1)).T * prob

    # multiply the y coord of the spots by the cell prob
    b = np.tile(xyz_spots[:, 1], (n, 1)).T * prob

    # multiply the z coord of the spots by the cell prob
    c = np.tile(xyz_spots[:, 2], (n, 1)).T * prob

    # aggregated x and y coordinate
    idx = spots.parent_cell_id
    x_agg = npg.aggregate(idx.ravel(), a.ravel(), size=len(N_c))
    y_agg = npg.aggregate(idx.ravel(), b.ravel(), size=len(N_c))
    z_agg = npg.aggregate(idx.ravel(), c.ravel(), size=len(N_c))

    # get the estimated cell centers
    x_bar = np.nan * np.ones(N_c.shape)
    y_bar = np.nan * np.ones(N_c.shape)
    z_bar = np.nan * np.ones(N_c.shape)

    x_bar[N_c > 0] = x_agg[N_c > 0] / N_c[N_c > 0]
    y_bar[N_c > 0] = y_agg[N_c > 0] / N_c[N_c > 0]
    z_bar[N_c > 0] = z_agg[N_c > 0] / N_c[N_c > 0]

    # cells with N_c = 0 will end up with x_bar = y_bar = np.nan
    xyz_bar_fitted = np.array(list(zip(x_bar.T, y_bar.T, z_bar.T)))

    # if you have a value for the estimated centroid use that, otherwise
    # use the initial (starting values) centroids
    ini_cent = cells.ini_centroids()
    xyz_bar = np.array(tuple(zip(*[ini_cent['x'], ini_cent['y'], ini_cent['z']])))

    # # sanity check. NaNs or Infs should appear together
    # assert np.all(np.isfinite(x_bar) == np.isfinite(y_bar))
    # use the fitted centroids where possible otherwise use the initial ones
    xyz_bar[np.isfinite(x_bar)] = xyz_bar_fitted[np.isfinite(x_bar)]
    return pd.DataFrame(xyz_bar, columns=['x', 'y', 'z'], dtype=np.float32)

