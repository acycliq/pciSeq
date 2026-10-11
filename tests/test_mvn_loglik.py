"""
Checks spots.mvn_loglik, the position term of the spot score.

It is the expected log density of a spot's position under the NIW posterior
of its cell (the-model/gaussian.md, section The expected log density), not a
plain gaussian log pdf. Two checks:

- against the closed form written out by hand here, with scipy's digamma, so
  the vectorised version and the formula on the docs page agree;
- in the limit kappa, nu -> infinity it must turn into scipy's
  multivariate_normal.logpdf with covariance V / nu, which is an independent
  implementation.
"""
import types

import numpy as np
import pytest
from scipy.special import digamma
from scipy.stats import multivariate_normal

from pciSeq.src.core.datatypes.spots import Spots


def _random_spd(rng, n, d):
    """n random symmetric positive-definite d-by-d matrices."""
    a = rng.standard_normal((n, d, d))
    return a @ np.transpose(a, (0, 2, 1)) + d * np.eye(d)


def _fake_cells(centroids, kappa, nu, V):
    """The attributes mvn_loglik reads, nothing else. centroids are always x, y, z."""
    vals, vecs = np.linalg.eigh(V)
    return types.SimpleNamespace(
        centroid=types.SimpleNamespace(values=centroids),
        kappa=kappa, nu=nu, scale=V, eig_vals=vals, eig_vecs=vecs,
    )


def _expected_log_density_by_hand(x, m, kappa, nu, V):
    """One spot, one cell, straight from the docs page."""
    d = len(x)
    Vinv = np.linalg.inv(V)
    E_logdet = sum(digamma((nu + 1 - i) / 2) for i in range(1, d + 1)) + d * np.log(2) - np.linalg.slogdet(V)[1]
    E_maha = d / kappa + nu * (x - m) @ Vinv @ (x - m)
    return -d / 2 * np.log(2 * np.pi) + 0.5 * E_logdet - 0.5 * E_maha


@pytest.mark.parametrize("d", [2, 3])
def test_matches_the_closed_form(d):
    rng = np.random.default_rng(0)
    nC, nS = 12, 80
    centroids = np.zeros((nC, 3))
    centroids[:, :d] = rng.standard_normal((nC, d))
    kappa = rng.uniform(10, 200, nC)
    nu = rng.uniform(10, 200, nC)
    V = _random_spd(rng, nC, d) * 5
    cells = _fake_cells(centroids, kappa, nu, V)

    cell_label = rng.integers(0, nC, size=nS)
    data = np.zeros((nS, 3))
    data[:, :d] = rng.standard_normal((nS, d))

    spots = Spots.__new__(Spots)   # no __init__: mvn_loglik only reads its arguments
    got = spots.mvn_loglik(data, cell_label, cells, is3D=(d == 3))

    ref = np.array([
        _expected_log_density_by_hand(data[i, :d], centroids[c, :d], kappa[c], nu[c], V[c])
        for i, c in enumerate(cell_label)
    ])
    np.testing.assert_allclose(got, ref, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("d", [2, 3])
def test_turns_into_the_gaussian_logpdf_when_the_posterior_is_certain(d):
    """kappa, nu -> inf: E[log N] -> log N(x | m, V / nu). Checked against scipy."""
    rng = np.random.default_rng(1)
    nC, nS = 8, 50
    big = 1e7
    centroids = np.zeros((nC, 3))
    centroids[:, :d] = rng.standard_normal((nC, d))
    covs = _random_spd(rng, nC, d)
    kappa = np.full(nC, big)
    nu = np.full(nC, big)
    V = covs * big
    cells = _fake_cells(centroids, kappa, nu, V)

    cell_label = rng.integers(0, nC, size=nS)
    data = np.zeros((nS, 3))
    data[:, :d] = rng.standard_normal((nS, d))

    spots = Spots.__new__(Spots)
    got = spots.mvn_loglik(data, cell_label, cells, is3D=(d == 3))

    ref = np.array([
        multivariate_normal(centroids[c, :d], covs[c]).logpdf(data[i, :d])
        for i, c in enumerate(cell_label)
    ])
    # the gap closes like d^2 / nu, so 1e-6 is plenty at nu = 1e7
    np.testing.assert_allclose(got, ref, rtol=1e-6, atol=1e-6)
