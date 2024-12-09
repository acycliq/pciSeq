"""Statistical calculation utilities."""
import numpy as np
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
