"""Geometric calculation utilities for cell and spot analysis."""
import numpy as np
import pandas as pd
from typing import List, Tuple, Union, Optional
from numpy.typing import NDArray
from scipy.sparse import coo_matrix
import logging

# Configure logging
geometry_logger = logging.getLogger(__name__)


def gaussian_contour(mu: Union[List, np.ndarray],
                     cov: np.ndarray,
                     sdwidth: float = 3) -> np.ndarray:
    """Draw an ellipsoid for a given covariance matrix and mean vector.

    Args:
        mu: Mean vector [x, y]
        cov: 2x2 covariance matrix
        sdwidth: Number of standard deviations to include (default: 3)

    Returns:
        Array of contour points (N x 2)

    Example:
        cov_1 = [[1, 0.5], [0.5, 1]]
        means_1 = [1, 1]
        ellipsis_1 = gaussian_contour(means_1, cov_1)
        plt.plot(ellipsis_1[0], ellipsis_1[1])
    """
    mu = np.array(mu)

    # Generate points around unit circle
    npts = 40
    tt = np.linspace(0, 2 * np.pi, npts)
    ap = np.zeros((2, npts))
    ap[0, :] = np.cos(tt)
    ap[1, :] = np.sin(tt)

    # Transform by eigenvalues/vectors of covariance
    eigvals, eigvecs = np.linalg.eig(cov)
    eigvals = sdwidth * np.sqrt(eigvals)
    eigvals = eigvals * np.eye(2)

    # Generate ellipse points
    vd = eigvecs.dot(eigvals)
    out = vd.dot(ap) + mu.reshape(2, -1)

    return np.array(list(zip(*out)), dtype=np.float32)


def gaussian_ellipsoid_props_OLD(cov: np.ndarray,
                             sdwidth: float = 3) -> Tuple[List[float], List[float]]:
    """Get the scaling and rotation of an ellipsoid.

    Args:
        cov: Covariance matrix
        sdwidth: Number of standard deviations (default: 3)

    Returns:
        Tuple containing:
            - List of scaling factors
            - List of rotation angles [theta_x, theta_y, theta_z]
    """
    tol = 1.0e-10
    cov = np.where(cov < tol, 0, cov)
    eigvals, eigvecs = np.linalg.eig(cov)
    scaling = sdwidth * np.sqrt(eigvals)
    rotation = euler_angles(eigvecs.T)
    return scaling.tolist(), rotation


def euler_angles_OLD(r: np.ndarray) -> List[float]:
    """Calculate Euler angles from rotation matrix.

    Args:
        r: 3x3 rotation matrix

    Returns:
        List of angles [theta_x, theta_y, theta_z]
    """
    theta_x = np.arctan2(r[2, 1], r[2, 2])
    theta_y = np.arcsin(r[2, 0])
    theta_z = np.arctan2(r[1, 0], r[0, 0])
    return [theta_x, theta_y, theta_z]


def gaussian_ellipsoid_props(cov: np.ndarray,
                             sdwidth: float = 3) -> Tuple[List[float], List[float]]:
    # Perform Eigen decomposition (PCA)
    tol = 1.0e-10
    cov = np.where(cov < tol, 0, cov)
    eigvals, eigvecs = np.linalg.eig(cov)

    # Sort eigenvalues and eigenvectors in descending order
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]  # Sort the eigenvalues
    eigvecs = eigvecs[:, idx]  # Sort the eigenvectors

    scaling = sdwidth * np.sqrt(eigvals)

    # Ensure eigenvectors are orthonormal (should already be, but this fixes small numerical issues)
    eigvecs, _ = np.linalg.qr(eigvecs)  # QR decomposition to make sure they are orthogonal

    # Ensure a right-handed coordinate system
    if np.linalg.det(eigvecs) < 0:
        eigvecs[:, -1] *= -1  # Flip the last eigenvector

    # Extract Euler angles from the rotation matrix
    r = eigvecs.T
    rotation = euler_angles(r)
    return scaling.tolist(), rotation


def euler_angles(r: NDArray[np.float32]) -> List[float]:
    """Extract Euler angles from a 3D rotation matrix using ZYX convention.

    Decomposes a rotation matrix into three sequential rotations around fixed
    axes (intrinsic rotations) following the ZYX (yaw-pitch-roll) convention:
    1. Rotation around Z axis (yaw/heading) - ψ (psi)
    2. Rotation around Y axis (pitch/attitude) - θ (theta)
    3. Rotation around X axis (roll/bank) - φ (phi)

    Parameters
    ----------
    r : NDArray[float32], shape (3, 3)
        Rotation matrix R where:
            R = Rz(ψ) @ Ry(θ) @ Rx(φ)
        Must be orthogonal with determinant 1

    Returns
    -------
    tuple of float
        Euler angles in radians (φ, θ, ψ) where:
            - φ (phi): rotation around X axis [-π/2, π/2]
            - θ (theta): rotation around Y axis [-π/2, π/2]
            - ψ (psi): rotation around Z axis [-π, π]
    """
    # Input validation
    if r.shape != (3, 3):
        raise ValueError("Input matrix must be 3x3")
    if not np.allclose(r @ r.T, np.eye(3), atol=1e-6):
        raise ValueError("Input matrix must be orthogonal")
    if not np.isclose(np.linalg.det(r), 1, atol=1e-6):
        raise ValueError("Input matrix must have determinant 1")

    # Extract Euler angles
    # Note: np.arctan2 returns angles in [-π, π]
    phi = np.arctan2(r[2, 1], r[2, 2])  # X rotation (roll)
    theta = np.arcsin(-r[2, 0])  # Y rotation (pitch)
    psi = np.arctan2(r[1, 0], r[0, 0])  # Z rotation (yaw)

    # Check for gimbal lock.
    # Gimbal lock occurs when θ = ±π/2, making φ and ψ indistinguishable.
    # if np.isclose(abs(theta), np.pi / 2, atol=1e-6):
    #     geometry_logger.warning("Warning: Gimbal lock detected (pitch = ±90°)")

    return [phi, theta, psi]


def adjust_for_anisotropy(
        spots: pd.DataFrame,
        voxel_size: Tuple[float, float, float]
) -> pd.DataFrame:
    """Adjust spot coordinates for non-cubic voxel dimensions.

    Transforms spatial coordinates to account for anisotropic voxels while
    preserving gene names and z-plane information. Applies scaling based on
    physical voxel dimensions.

    Parameters
    ----------
    spots : pd.DataFrame
        Spot data with columns:
            - gene_name : str, gene identifier
            - x : float, x coordinate
            - y : float, y coordinate
            - z_plane : float, z-plane index

    voxel_size : tuple of float, length 3
        Physical voxel dimensions (dx, dy, dz) in consistent units.
        Example: (0.5, 0.5, 2.0) for 0.5μm x 0.5μm x 2.0μm voxels

    Returns
    -------
    pd.DataFrame
        Adjusted spot data with columns:
            - gene_name : str
            - x : float32, scaled x coordinate
            - y : float32, scaled y coordinate
            - z : float32, scaled z coordinate
            - z_plane : float32, original z-plane index
    """

    # Extract columns for transformation
    gene_col = spots.gene_name.values[:, None]  # Add dimension for hstack
    z_plane = spots.z_plane.values[:, None].astype(np.float32)
    score = spots.score.values[:, None].astype(np.float32)
    intensity = spots.intensity.values[:, None].astype(np.float32)

    # Get coordinates for scaling
    coords = spots[['x', 'y', 'z_plane']].values

    # Apply anisotropic scaling
    scaled_coords = anisotropy_calc(coords, voxel_size)

    # Combine all columns
    data_adj = np.hstack([
        gene_col,  # Gene names
        scaled_coords,  # Scaled coordinates
        z_plane,  # Original z-plane indices
        score,
        intensity
    ])

    # Create DataFrame with correct column names and data types.
    # Note for future: consider preserving the original index from the input 'spots' DataFrame.
    # Keeping the index can be useful for filtering spots before calling pciSeq.fit().
    # The line below creates a new index, which can make it harder to exclude spots
    # if they were removed or cleaned earlier in the codebase.
    # Should be safe, but worth testing.
    return pd.DataFrame(
        data=data_adj,
        columns=['gene_name', 'x', 'y', 'z', 'z_plane', 'score', 'intensity']
    ).astype({
        'gene_name': str,
        'x': np.float32,
        'y': np.float32,
        'z': np.float32,
        'z_plane': np.float32,
        'score': np.float32,
        'intensity': np.float32,
    })


def anisotropy_calc(data: np.ndarray,
                    voxel_size: Tuple[float, float, float],
                    inverse = False) -> np.ndarray:
    """Calculate anisotropic scaling for spot coordinates.

    Adjusts coordinates to account for different voxel dimensions in x, y, and z.
    All dimensions are scaled relative to the x dimension.

    Args:
        data: Array of shape (N, 3) containing spot coordinates [x, y, z]
        voxel_size: Physical dimensions of voxels as (x, y, z) in same units
        inverse: If True, return the original coordinates. If False, return the scaled coordinates. Default is False.

    Returns:
        np.ndarray: Scaled coordinates of shape (N, 3)

    Notes:
        - X dimension is used as reference (Sx = 1)
        - Y and Z are scaled relative to X
        - Output maintains same shape as input
        - Uses float32 for memory efficiency
    """
    # Unpack voxel dimensions
    x, y, z = voxel_size

    # Calculate scaling factors relative to x dimension
    Sx = x / x  # Always 1
    Sy = y / x  # Scale y relative to x
    Sz = z / x  # Scale z relative to x

    # Create anisotropic scaling matrix
    scaling_matrix = np.array([
        [Sx, 0, 0],
        [0, Sy, 0],
        [0, 0, Sz]
    ], dtype=np.float32)

    if not inverse:
        # Apply scaling: matrix multiplication then transpose back
        # data.T converts from (N,3) to (3,N) for matrix multiplication
        # final .T converts back to (N,3)
        return scaling_matrix.dot(data.T).T
    elif inverse:
        # does the inverse of the above.
        # Coordinates are already scaled relative to x dimension, and
        # we want to return the original coordinates.
        # In most cases where the voxel size differs only on z (ie Sx=Sy but Sz!=Sx)
        # that will be equivalent as dividing z by Sz/Sx
        S = scaling_matrix
        W = np.linalg.inv(S.T.dot(S)).dot(S.T)
        return W.dot(data.T).T
    else:
        raise ValueError('inverse must be True or False')





def get_img_shape(coo: List[coo_matrix]) -> List[int]:
    """Get image dimensions from list of sparse matrices.

    Args:
        coo: List of sparse coordinate matrices

    Returns:
        List of dimensions [n_planes, height, width]

    Raises:
        AssertionError: If matrices have different shapes
    """
    n = len(coo)
    img_shape = set([d.shape for d in coo])
    assert len(img_shape) == 1, 'pages do not have the same shape'
    img_shape = img_shape.pop()
    return [n, img_shape[0], img_shape[1]]
