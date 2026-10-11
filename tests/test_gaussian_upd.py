"""
gaussian_upd, the cell geometry step: the NIW posterior over each cell's centre
and covariance, the-model/gaussian.md.

Checked against a plain python loop over spots that does the update one cell at
a time, so the vectorised numpy_groupies version is compared with something
written straight from the docs page and nothing else.
"""
import numpy as np
import pandas as pd
import pytest

from pciSeq.src.core.main import VarBayes
from pciSeq.src.core.utils.elbo import geometry_kl, calc_elbo


def _one_pass(vb):
    vb.initialise_state()
    for step in (vb.geneCount_upd, vb.rho_upd, vb.eta_upd, vb.theta_upd,
                 vb.gamma_upd, vb.cell_to_cellType, vb.spots_to_cell_numba):
        step()


def _niw_by_hand(vb, c):
    """kappa_hat, m_hat, nu_hat, V_hat of cell c from a python loop over spots."""
    cfg = vb.config
    kappa0 = cfg['cell_centroid_prior']['default']
    nu0 = cfg['cell_cov_prior']['default']
    d = vb.cells.dim
    r = vb.cells.mcr
    m0 = vb.cells.ini_centroids().values[c, :d].astype(np.float64)

    xyz = vb.spots.xyz_coords[:, :d].astype(np.float64)
    ids = vb.spots.parent_cell_id[:, :-1]
    prob = vb.spots.parent_cell_prob[:, :-1]

    n = 0.0
    sx = np.zeros(d)
    for s in range(xyz.shape[0]):
        for j in range(ids.shape[1]):
            if ids[s, j] == c:
                n += prob[s, j]
                sx += prob[s, j] * xyz[s]
    x_bar = sx / n if n > 0 else m0.copy()

    S = np.zeros((d, d))
    for s in range(xyz.shape[0]):
        for j in range(ids.shape[1]):
            if ids[s, j] == c:
                dev = xyz[s] - x_bar
                S += prob[s, j] * np.outer(dev, dev)

    kappa_hat = kappa0 + n
    m_hat = (kappa0 * m0 + n * x_bar) / kappa_hat
    nu_hat = nu0 + n
    V0 = (nu0 - d - 1) * r * r * np.eye(d)
    V_hat = V0 + S + kappa0 * n / (kappa0 + n) * np.outer(x_bar - m0, x_bar - m0)
    return kappa_hat, m_hat, nu_hat, V_hat


def test_runs_and_moves_the_cells(minimal_varbayes):
    vb = minimal_varbayes
    _one_pass(vb)
    cov_before = vb.cells.cov.copy()
    centroid_before = vb.cells.centroid.values.copy()

    vb.gaussian_upd()

    assert np.isfinite(vb.cells.cov).all()
    assert np.isfinite(vb.cells.centroid.values[1:]).all()
    assert vb.cells.cov.shape == cov_before.shape
    assert vb.cells.centroid.values.shape == centroid_before.shape
    assert not np.allclose(vb.cells.cov, cov_before), 'cov did not move at all'
    assert not np.allclose(vb.cells.centroid.values[1:], centroid_before[1:]), 'centroids did not move'


def test_matches_the_update_done_by_hand(minimal_varbayes):
    vb = minimal_varbayes
    _one_pass(vb)
    vb.gaussian_upd()
    d = vb.cells.dim
    for c in range(1, vb.nC):
        kappa_hat, m_hat, nu_hat, V_hat = _niw_by_hand(vb, c)
        np.testing.assert_allclose(vb.cells.kappa[c], kappa_hat, rtol=1e-5)
        np.testing.assert_allclose(vb.cells.centroid.values[c, :d], m_hat, rtol=1e-5)
        np.testing.assert_allclose(vb.cells.nu[c], nu_hat, rtol=1e-5)
        np.testing.assert_allclose(vb.cells.scale[c], V_hat, rtol=1e-4)
        np.testing.assert_allclose(vb.cells.cov[c], V_hat / (nu_hat - d - 1), rtol=1e-4)


def test_the_background_stays_on_the_prior(minimal_varbayes):
    """Row 0 is not a cell. It gets no spots and its posterior is its prior."""
    vb = minimal_varbayes
    _one_pass(vb)
    vb.gaussian_upd()
    cfg = vb.config
    assert vb.cells.kappa[0] == cfg['cell_centroid_prior']['default']
    assert vb.cells.nu[0] == cfg['cell_cov_prior']['default']
    np.testing.assert_allclose(vb.cells.scale[0], vb.cells.ini_scale()[0])
    np.testing.assert_allclose(vb.cells.cov[0], vb.cells.ini_cov()[0], rtol=1e-6)


def test_a_cell_with_no_spots_stays_on_the_prior(minimal_varbayes):
    vb = minimal_varbayes
    _one_pass(vb)
    c = 3
    prob = vb.spots.parent_cell_prob.copy()
    prob[vb.spots.parent_cell_id == c] = 0.0
    vb.spots.parent_cell_prob = prob
    vb.gaussian_upd()
    cfg = vb.config
    d = vb.cells.dim
    assert vb.cells.kappa[c] == cfg['cell_centroid_prior']['default']
    assert vb.cells.nu[c] == cfg['cell_cov_prior']['default']
    np.testing.assert_allclose(vb.cells.centroid.values[c, :d], vb.cells.ini_centroids().values[c, :d])
    np.testing.assert_allclose(vb.cells.scale[c], vb.cells.ini_scale()[c])


def test_the_prior_scale_puts_the_prior_mean_at_r_squared(minimal_varbayes):
    """V0 = (nu0 - d - 1) r^2 I so that E[Sigma] = V0 / (nu0 - d - 1) = r^2 I (Zhang 2021)."""
    vb = minimal_varbayes
    nu0 = vb.config['cell_cov_prior']['default']
    d = vb.cells.dim
    np.testing.assert_allclose(vb.cells.ini_scale() / (nu0 - d - 1), vb.cells.ini_cov(), rtol=1e-6)


def test_geometry_kl_is_zero_on_the_prior_and_positive_after(minimal_varbayes):
    vb = minimal_varbayes
    _one_pass(vb)
    assert geometry_kl(vb) == pytest.approx(0.0, abs=1e-6)
    vb.gaussian_upd()
    assert geometry_kl(vb) > 0
    assert np.isfinite(calc_elbo(vb))


def test_cell_cov_prior_must_exceed_d_plus_one(rng, base_opts, minimal_varbayes):
    """nu0 <= d + 1 has no prior mean, so it is refused up front."""
    vb = minimal_varbayes
    bad = {**vb.config, 'cell_cov_prior': {'default': 4}}
    cells_df = pd.DataFrame({'label': [1, 2], 'area': [100, 100], 'x0': [1, 2], 'y0': [1, 2], 'z0': [1, 2]})
    from pciSeq.src.core.datatypes.cells import Cells
    with pytest.raises(ValueError):
        Cells(cells_df, bad)


def test_two_dimensional_run(rng, base_opts):
    """2D used to crash in the old scatter_matrix (a 3x3 write into a 2x2). Same
    synthetic data as the fixture, flattened, and checked by hand as well."""
    gene_names = [f"Gene_{i}" for i in range(10)]
    scref_df = pd.DataFrame(rng.random((10, 3)) * 100, index=gene_names, columns=["A", "B", "C"])
    scref_df.index.name = "gene_name"
    scref_df.columns.name = "class_name"
    cells_df = pd.DataFrame({
        "label": [1, 2, 3, 4, 5, 6, 7, 8],
        "area": [100, 120, 90, 110, 95, 105, 115, 100],
        "x0": [10, 20, 30, 40, 10, 20, 30, 40],
        "y0": [10, 20, 30, 40, 50, 50, 50, 50],
        "z0": [0] * 8,
    })
    spots_df = pd.DataFrame({
        "x": rng.uniform(0, 50, 30),
        "y": rng.uniform(0, 50, 30),
        "z": np.zeros(30),
        "plane_id": np.zeros(30, dtype=int),
        "label": rng.choice([1, 2, 3, 4, 5, 6, 7, 8], 30),
        "gene_name": rng.choice(gene_names, 30),
        "score": rng.uniform(0.5, 1.0, 30),
        "intensity": rng.uniform(100, 1000, 30),
    })
    opts = {**base_opts, "is3D": False, "img_dim": {"w": 100, "h": 100, "n_planes": 1},
            "label_map": {}, "MisreadDensity": {"default": 1e-6}, "rRho": 1000.0,
            "rTheta": 25.0, "mrf_beta": 1.0}
    vb = VarBayes(spots_df=spots_df, cells_df=cells_df, scRNAseq=scref_df, config=opts)
    assert vb.cells.dim == 2
    _one_pass(vb)
    vb.gaussian_upd()
    assert vb.cells.scale.shape == (vb.nC, 2, 2)
    assert np.isfinite(vb.cells.cov).all()
    for c in range(1, vb.nC):
        kappa_hat, m_hat, nu_hat, V_hat = _niw_by_hand(vb, c)
        np.testing.assert_allclose(vb.cells.centroid.values[c, :2], m_hat, rtol=1e-5)
        np.testing.assert_allclose(vb.cells.scale[c], V_hat, rtol=1e-4)
    # and the spot scores still work on the 2D geometry
    vb.spots_to_cell_numba()
    assert np.isfinite(vb.spots.mvn_loglik_arr[:, :-1]).all()


def test_cellData_keeps_the_segmentation_centroid(minimal_varbayes):
    """X, Y, Z is the fitted centre, X0, Y0, Z0 where the segmentation put the
    cell. After one geometry update the two differ for some cell."""
    from pciSeq.src.core.summary import collect_data
    vb = minimal_varbayes
    _one_pass(vb)
    vb.gaussian_upd()

    cellData, _ = collect_data(vb.cells, vb.spots, vb.genes, vb.config['is3D'])
    assert list(cellData.columns)[:7] == ['Cell_Num', 'X', 'Y', 'Z', 'X0', 'Y0', 'Z0']

    m0 = vb.cells.ini_centroids().loc[cellData['Cell_Num'].values]
    m_hat = vb.cells.centroid.loc[cellData['Cell_Num'].values]
    for col, src in (('X0', 'x'), ('Y0', 'y'), ('Z0', 'z')):
        np.testing.assert_allclose(cellData[col].values, m0[src].values, atol=1e-3)
    for col, src in (('X', 'x'), ('Y', 'y'), ('Z', 'z')):
        np.testing.assert_allclose(cellData[col].values, m_hat[src].values, atol=1e-3)
    assert not np.allclose(cellData[['X', 'Y']].values, cellData[['X0', 'Y0']].values), \
        'the fitted centre never left the segmentation centroid'
