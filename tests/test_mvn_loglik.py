"""
Checks mvn_loglik / multiple_logpdfs against scipy's standard multivariate_normal
as an independent reference.

`multiple_logpdfs` (and its caller `mvn_loglik`) evaluate the MVN log-pdf for many
(point, mean, covariance) triples at once, using the eigendecomposition of each
covariance. These tests check that result against scipy's
`multivariate_normal.logpdf`, one triple at a time.

The point is to check against an independent implementation: when we later
reimplement the maths (numba), we verify the new version against scipy here, not
just against the old code.
"""
import types

import numpy as np
import pytest
from scipy.stats import multivariate_normal

from pciSeq.src.core.datatypes.spots import Spots


def _random_spd(rng, n, d):
    """n random symmetric positive-definite d-by-d matrices."""
    a = rng.standard_normal((n, d, d))
    return a @ np.transpose(a, (0, 2, 1)) + d * np.eye(d)


@pytest.mark.parametrize("d", [2, 3])
def test_multiple_logpdfs_matches_scipy(d):
    rng = np.random.default_rng(0)
    n = 64
    x = rng.standard_normal((n, d))
    means = rng.standard_normal((n, d))
    covs = _random_spd(rng, n, d)
    vals, vecs = np.linalg.eigh(covs)

    spots = Spots.__new__(Spots)  # no __init__: multiple_logpdfs uses only its args
    got = spots.multiple_logpdfs(x, means, covs, vals, vecs)

    ref = np.array([multivariate_normal(means[i], covs[i]).logpdf(x[i]) for i in range(n)])
    np.testing.assert_allclose(got, ref, rtol=1e-10, atol=1e-10)


def test_mvn_loglik_matches_scipy_3d():
    """mvn_loglik gathers each spot's cell parameters, then calls multiple_logpdfs."""
    rng = np.random.default_rng(1)
    nC, nS, d = 12, 80, 3
    centroids = rng.standard_normal((nC, d))
    covs = _random_spd(rng, nC, d)
    vals, vecs = np.linalg.eigh(covs)

    # minimal stand-in for the Cells object (only the attributes mvn_loglik reads)
    cells = types.SimpleNamespace(
        centroid=types.SimpleNamespace(values=centroids),
        cov=covs,
        eig_vals=vals,
        eig_vecs=vecs,
    )
    cell_label = rng.integers(0, nC, size=nS)
    data = rng.standard_normal((nS, d))

    spots = Spots.__new__(Spots)
    got = spots.mvn_loglik(data, cell_label, cells, is3D=True)

    ref = np.array([
        multivariate_normal(centroids[c], covs[c]).logpdf(data[i])
        for i, c in enumerate(cell_label)
    ])
    np.testing.assert_allclose(got, ref, rtol=1e-10, atol=1e-10)
