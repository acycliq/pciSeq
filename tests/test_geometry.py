"""gaussian_ellipsoid_props and ellipsoid_to_cov must be inverses.

sphere_scale and sphere_rotation are the only place cellData keeps the fitted
covariance, so the six numbers have to encode the whole 3x3 matrix.
"""
import numpy as np
import pytest

from pciSeq.src.core.utils.geometry import gaussian_ellipsoid_props, ellipsoid_to_cov


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
