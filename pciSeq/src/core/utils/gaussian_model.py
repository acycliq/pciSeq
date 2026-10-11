"""The gaussian each cell's spots are taken to be scattered around.

Its centre m_c and covariance Sigma_c get a Normal-Inverse-Wishart prior and
the posterior comes out NIW again. The maths is on the cell geometry page,
website/docs/the-model/gaussian.md, and the symbols here follow it:
n_bar, x_bar, S are the soft count, mean and scatter of a cell's spots,
kappa, nu, V the NIW parameters, hat for the posterior, 0 for the prior.
"""

import logging
import numpy as np
import numpy_groupies as npg
import opt_einsum as oe
from scipy.special import digamma, multigammaln

logger = logging.getLogger(__name__)


def soft_stats(spots, cells, d):
    """Soft count, soft mean and soft scatter of the spots around every cell.

        n_bar[c] = sum_s z_sc
        x_bar[c] = sum_s z_sc x_s / n_bar[c]
        S[c]     = sum_s z_sc (x_s - x_bar[c]) (x_s - x_bar[c])^T

    z_sc is parent_cell_prob. The background column is left out, so cell 0
    gets no spots and stays on its prior. A cell with n_bar = 0 gets x_bar
    equal to its segmentation centroid, which makes the centre shift term of
    the scale update zero for it. Sums run in float64: squared coordinates are
    big enough for float32 to lose the spread.

    Parameters
    ----------
    spots : Spots
    cells : Cells
    d : int
        2 or 3, the number of spatial dimensions. z is dropped in 2D.

    Returns
    -------
    n_bar : (nC,) float64
    x_bar : (nC, d) float64
    S : (nC, d, d) float64
    """
    nC = cells.nC
    ids = spots.parent_cell_id[:, :-1].ravel()
    w = spots.parent_cell_prob[:, :-1].ravel().astype(np.float64)
    nNb = spots.parent_cell_id.shape[1] - 1
    # one row per (spot, neighbour) pair, same order as the ravel above
    X = np.repeat(spots.xyz_coords[:, :d].astype(np.float64), nNb, axis=0)

    n_bar = npg.aggregate(ids, w, size=nC)

    x_bar = np.empty((nC, d))
    for j in range(d):
        x_bar[:, j] = npg.aggregate(ids, w * X[:, j], size=nC)
    has_spots = n_bar > 0
    x_bar[has_spots] /= n_bar[has_spots, None]
    m0 = cells.ini_centroids().values[:, :d].astype(np.float64)
    x_bar[~has_spots] = m0[~has_spots]

    # scatter about the soft mean, not about the raw second moment: the two
    # differ by n_bar x_bar x_bar^T which is huge next to the spread
    dev = X - x_bar[ids]
    S = np.empty((nC, d, d))
    for j in range(d):
        for k in range(j, d):
            S[:, j, k] = npg.aggregate(ids, w * dev[:, j] * dev[:, k], size=nC)
            S[:, k, j] = S[:, j, k]
    return n_bar, x_bar, S


def expected_logdet_precision(nu, eig_vals):
    """E[log|Lambda|] for Lambda ~ Wishart(V^-1, nu), one value per cell.

        sum_{i=1..d} psi((nu + 1 - i) / 2) + d log 2 - log|V|

    Parameters
    ----------
    nu : (nC,)
    eig_vals : (nC, d)
        eigenvalues of V, so log|V| is the sum of their logs
    """
    d = eig_vals.shape[1]
    i = np.arange(1, d + 1)
    psi_sum = digamma((nu[:, None] + 1 - i) / 2).sum(axis=1)
    return psi_sum + d * np.log(2) - np.log(eig_vals).sum(axis=1)


def niw_kl(m_hat, kappa_hat, V_hat, nu_hat, m0, kappa0, V0, nu0):
    """KL(q || p) between the NIW posterior and the NIW prior, one value per cell.

    The two pieces, the gaussian part averaged over Sigma and the inverse
    Wishart part, are written out on the cell geometry page under the ELBO.
    All inputs are per cell except kappa0 and nu0 which are scalars.
    """
    nC, d = m_hat.shape
    m_hat = m_hat.astype(np.float64)
    m0 = m0.astype(np.float64)
    V_hat = V_hat.astype(np.float64)
    V0 = V0.astype(np.float64)
    nu_hat = nu_hat.astype(np.float64)
    kappa_hat = kappa_hat.astype(np.float64)

    vals, vecs = np.linalg.eigh(V_hat)
    logdet_hat = np.log(vals).sum(axis=1)
    logdet0 = np.linalg.slogdet(V0)[1]
    V_hat_inv = oe.contract('cij,cj,ckj->cik', vecs, 1.0 / vals, vecs)

    dm = m_hat - m0
    quad = oe.contract('ci,cij,cj->c', dm, V_hat_inv, dm)
    kl_N = (d / 2 * np.log(kappa_hat / kappa0) - d / 2 + d * kappa0 / (2 * kappa_hat)
            + kappa0 * nu_hat / 2 * quad)

    E_logdet_L = expected_logdet_precision(nu_hat, vals)
    tr = oe.contract('cij,cji->c', V0, V_hat_inv)
    kl_IW = (nu_hat / 2 * logdet_hat - nu0 / 2 * logdet0
             - (nu_hat - nu0) * d / 2 * np.log(2)
             - multigammaln(nu_hat / 2, d) + multigammaln(nu0 / 2, d)
             + (nu_hat - nu0) / 2 * E_logdet_L
             + nu_hat / 2 * (tr - d))
    return kl_N + kl_IW
