"""gaussian_ellipsoid_props and ellipsoid_to_cov must be inverses.

sphere_scale and sphere_rotation are the only place cellData keeps the fitted
covariance, so the six numbers have to encode the whole 3x3 matrix.
"""
import numpy as np
import pytest

from pciSeq.src.core.utils.geometry import (gaussian_ellipsoid_props, ellipsoid_to_cov,
                                             plane_of, plane_to_z, anisotropy_calc)


@pytest.mark.parametrize('seed', range(20))
def test_scale_and_rotation_rebuild_the_covariance(seed):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(3, 3))
    cov = (a @ a.T * 20).astype(np.float32)
    scale, rotation = gaussian_ellipsoid_props(cov)
    assert not np.isnan(rotation).any()
    assert scale == sorted(scale, reverse=True)
    np.testing.assert_allclose(ellipsoid_to_cov(scale, rotation), cov, rtol=1e-4, atol=1e-3)


def test_negative_off_diagonals_are_kept():
    # used to be zeroed by the noise filter, which turned the matrix into another one
    cov = np.array([[9.0, -4.0, 0.0], [-4.0, 9.0, 0.0], [0.0, 0.0, 4.0]])
    scale, rotation = gaussian_ellipsoid_props(cov)
    np.testing.assert_allclose(sorted(scale, reverse=True), 3 * np.sqrt([13.0, 5.0, 4.0]))
    # the angles are rounded to 7 digits, hence the tolerance
    np.testing.assert_allclose(ellipsoid_to_cov(scale, rotation), cov, atol=1e-5)


def test_fixed_gaussian_is_a_sphere_of_three_mcr():
    mcr = 10.328247
    cov = mcr * mcr * np.eye(3, dtype=np.float32)
    scale, rotation = gaussian_ellipsoid_props(cov)
    np.testing.assert_allclose(scale, [3 * mcr] * 3, rtol=1e-6)
    assert not np.isnan(rotation).any()
    np.testing.assert_allclose(ellipsoid_to_cov(scale, rotation), cov, rtol=1e-5, atol=1e-6)


# The same numbers are in pciSeq_viewer tests/planes.test.js. If one side changes
# the rule the other side's test still holds the old numbers, so change both.
# (voxel_size, scaled z, plane)
PLANE_CASES = [
    ([0.28, 0.28, 0.7], 0.0, 0),
    ([0.28, 0.28, 0.7], 2.875, 1),     # spot 2293181 on espio_demo, plane_id 1
    ([0.28, 0.28, 0.7], 52.49, 20),    # just below plane 21
    ([0.28, 0.28, 0.7], 52.5, 21),     # exactly on plane 21
    ([0.28, 0.28, 0.7], 53.91, 21),    # plane 21.56, floor not round
    ([1, 1, 1], 0.9998, 0),            # outside the 1e-4 tolerance, still 0
    ([1, 1, 1], 0.99995, 1),           # inside it, taken as on plane 1
    ([1, 1, 1], 7.5, 7),
]


@pytest.mark.parametrize('voxel_size, z, plane', PLANE_CASES)
def test_plane_of_is_the_floor(voxel_size, z, plane):
    assert plane_of(z, voxel_size) == plane


@pytest.mark.parametrize('voxel_size', [[0.28, 0.28, 0.7], [0.147, 0.147, 0.9]])
def test_plane_to_z_and_back(voxel_size):
    # exactly on a plane, even when the ratio is not exact in floating point
    planes = np.arange(100)
    assert (plane_of(plane_to_z(planes, voxel_size), voxel_size) == planes).all()
    assert all(plane_of(plane_to_z(p, voxel_size), voxel_size) == p for p in range(100))


# fractional planes, as a spot caller gives them; the viewer test uses the same
ROUND_TRIP_PLANES = [0.0, 0.3, 1.15, 20.999, 21.0, 21.56, 88.85]


@pytest.mark.parametrize('voxel_size', [[0.28, 0.28, 0.7], [0.147, 0.147, 0.9], [1, 1, 4]])
def test_z_plane_to_z_and_back_before_the_floor(voxel_size):
    """The round trip without the floor, which cannot be undone, and then the floor
    on its own. plane_to_z has to give the z that anisotropy_calc, which really
    builds the spot z, gives; dividing back has to give the plane; and plane_of has
    to be the floor of it."""
    zp = np.array(ROUND_TRIP_PLANES)
    z = plane_to_z(zp, voxel_size)
    xyz = np.stack([np.zeros_like(zp), np.zeros_like(zp), zp], axis=1)
    np.testing.assert_allclose(z, anisotropy_calc(xyz, voxel_size)[:, 2], rtol=1e-6)
    vx, _, vz = voxel_size
    np.testing.assert_allclose(z * vx / vz, zp, rtol=1e-12, atol=1e-12)
    assert (plane_of(z, voxel_size) == np.floor(zp)).all()

