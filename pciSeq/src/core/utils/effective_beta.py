"""
Per-(cell, class) adaptive cap on the MRF strength used in the cell-typing
step. See pciSeq_model/mrf_inflection_point.tex for the full derivation.

For each cell c and each real class k, computes
    beta*_{c,k} = ( |D_k|_struct(c) - bonus_k(c) - log(pi_k / pi_Zero) ) / |N_c|
with
    |D_k|_struct = r_gamma * sum_g log[ (r_gamma + mu*xi) / (r_gamma + sigma*xi) ]
    bonus_k      = sum_g N_{c,g} * log[ mu/sigma * (r_gamma + sigma*xi) / (r_gamma + mu*xi) ]
    xi_g         = A_c * eta_bar_g * theta_bar_{c,k}   (A_c approx 1)
    (mu = mean_expression_adj, which already carries the SpotReg regulariser
     sigma baked in; the Zero floor uses sigma directly = the Zero column of mu.)

The effective beta used per (c, k) is then
    min(config['mrf_beta'], max(0, beta*_{c,k}))
with the Zero column left at config['mrf_beta'] unchanged. Cells where the
data and baseline prior together already prefer class k over Zero get a
negative beta*_{c,k} and pass through the full configured beta.

Two implementations are provided, kept side-by-side so they can be
cross-checked: a plain numpy per-class loop and a numba-parallel fused
kernel that runs about 4x faster on a 24k-cell run.
"""

import numba as nb
import numpy as np


@nb.njit(parallel=True, fastmath=True, cache=True)
def _calc_beta_cap_numba_kernel(mu_real, eta_bar, theta_bar, theta_zero, N, log_pi_ratio,
                                r, sigma, n_nbr):
    """
    Fused triple-nested kernel for the per-(cell, class) MRF cap. Equivalent
    to ``compute_effective_beta`` but does the elementwise log / multiply /
    sum in a single streaming pass over (c, g, k) with no temporary arrays,
    and parallelises across cells.

    log_pi_ratio is the (K_real,) vector log(pi_k / pi_Zero). Under a uniform
    baseline prior it is identically zero and contributes nothing; under a
    non-uniform prior (cell_type_weights set, or 'weighted' mode where the
    Dirichlet posterior evolves) it shifts each class's cap by a constant.

    Returns the raw beta_cap matrix of shape (nC, K_real); the np.where cap
    is applied by the caller.
    """
    nC, G = N.shape
    K_real = mu_real.shape[1]

    # Pre-bake log(mu/sigma): mu already carries the SpotReg regulariser (baked
    # into mean_expression_adj), so no extra +sigma on the class mean here.
    log_mu_ratio = np.empty((G, K_real), dtype=np.float32)
    for g in range(G):
        for k in range(K_real):
            log_mu_ratio[g, k] = np.log(mu_real[g, k] / sigma)

    beta_cap = np.empty((nC, K_real), dtype=np.float32)
    for c in nb.prange(nC):
        tz_c = theta_zero[c]                                  # theta_{c,Zero}
        for k in range(K_real):
            acc_log_ratio = np.float32(0.0)
            acc_bonus = np.float32(0.0)
            tb_ck = theta_bar[c, k]                           # theta_{c,k}
            log_theta_ratio = np.log(tb_ck / tz_c)           # log(theta_{c,k}/theta_{c,Zero})
            for g in range(G):
                num = r + mu_real[g, k] * (eta_bar[g] * tb_ck)             # r + m_k  (mu has eps baked in)
                den = r + sigma * (eta_bar[g] * tz_c)                       # r + m_0  (theta_{c,Zero})
                lr = np.log(num / den)
                acc_log_ratio += lr
                # bonus_g = N * log[ mu/sigma * theta_{c,k}/theta_{c,Zero} * den/num ]
                acc_bonus += N[c, g] * (log_mu_ratio[g, k] + log_theta_ratio - lr)
            beta_cap[c, k] = (r * acc_log_ratio - acc_bonus - log_pi_ratio[k]) / n_nbr
    return beta_cap

def beta_cap(obj):
    eta_bar = obj.genes.eta_bar
    theta_bar = obj.cells.theta_bar
    mu = obj.single_cell.mean_expression_adj.values
    spotReg = obj.config["SpotReg"]
    r = obj.config["rSpot"]
    beta_cfg = obj.config["mrf_beta"]

    a = r + np.einsum("gk,g,ck->cgk", mu[:,:-1], eta_bar, theta_bar[:,:-1])
    b = r + np.einsum("g,c->cg", spotReg * eta_bar, theta_bar[:,-1])   # (nC, nG)  Zero floor, theta_{c,Zero}
    theta_ratio = theta_bar[:, :-1] / theta_bar[:, -1][:, None]        # (nC, nK-1)  theta_{c,k}/theta_{c,Zero}

    struct = r * np.log(a / b[:, :, None]).sum(axis=1)

    bonus = (
        obj.cells.geneCount[:, :, None]
        * np.log((mu[:, :-1]) / spotReg * theta_ratio[:, None, :] * (b[:, :, None] / a))
    ).sum(axis=1)

    log_prior = np.log(obj.cellTypes.prior[:-1] / obj.cellTypes.prior[-1])

    beta_cap = (struct - bonus - log_prior) / obj.config["nNeighbors"]
    beta_cap = np.where(beta_cap < 0, beta_cfg, np.minimum(beta_cfg, beta_cap))

    out = np.full((obj.nC, obj.nK), beta_cfg, dtype=np.float32)
    out[:, :-1] = beta_cap

    return out


def compute_effective_beta(obj) -> np.ndarray:
    """
    Numpy reference implementation of the per-(cell, class) MRF cap.

    Reads the current variational state off ``obj`` (a VarBayes instance) and
    returns an (nC, nK) float32 matrix. The Zero column is left at the
    configured beta; the real-class columns are the cap.
    """
    r = obj.config['rSpot']
    spotReg = obj.config['SpotReg']
    n_nbr = obj.config['nNeighbors']
    beta_cfg = obj.config['mrf_beta']

    mu = obj.single_cell.mean_expression_adj.values       # (nG, nK)
    mu_real = mu[:, :-1]                                  # (nG, nK-1)
    eta_bar = np.asarray(obj.genes.eta_bar)               # (nG,)
    theta_bar_full = np.asarray(obj.cells.theta_bar)      # (nC, nK)
    theta_bar = theta_bar_full[:, :-1]                    # (nC, nK-1)  capture under each real class
    theta_zero = theta_bar_full[:, -1]                    # (nC,)       capture under the Zero class
    geneCount = np.asarray(obj.cells.geneCount)           # (nC, nG)

    # Baseline class prior correction log(pi_k / pi_Zero), refetched every
    # iteration so the formula stays correct under 'weighted' mode where
    # the Dirichlet posterior evolves. Identically zero under a uniform
    # prior with no cell_type_weights overrides.
    log_prior = np.asarray(obj.cellTypes.log_prior)       # (nK,)
    log_pi_ratio = log_prior[:-1] - log_prior[-1]         # (nK-1,)

    # The Zero floor uses theta_{c,Zero} (capture under the Zero class), which is
    # the same for every k, so build it once. The class-k mean uses theta_{c,k}.
    eta32 = eta_bar.astype(np.float32)
    tz = theta_zero.astype(np.float32)                    # theta_{c,Zero}
    xi_z = eta32[None, :] * tz[:, None]                   # (nC, nG)  eta * theta_{c,Zero}
    den = r + spotReg * xi_z                              # r + m_0,  Zero floor
    log_tz = np.log(tz)

    # per-class loop keeps peak memory at O(nC * nG) instead of O(nC * nG * nK)
    beta_cap = np.empty((obj.nC, obj.nK - 1), dtype=np.float32)
    for k in range(obj.nK - 1):
        mu_k = mu_real[:, k].astype(np.float32)
        tb = theta_bar[:, k].astype(np.float32)           # theta_{c,k}
        xi_k = eta32[None, :] * tb[:, None]               # (nC, nG)  eta * theta_{c,k}
        num = r + mu_k[None, :] * xi_k                    # r + m_k  (mu has eps baked in)
        D_struct = r * np.log(num / den).sum(axis=1)
        # m_k/m_0 = mu/eps * theta_{c,k}/theta_{c,Zero}  (eta cancels, theta does not)
        log_mk_over_m0 = np.log(mu_k / spotReg)[None, :] + (np.log(tb) - log_tz)[:, None]
        log_bonus_per_gene = log_mk_over_m0 + np.log(den / num)
        bonus = (geneCount * log_bonus_per_gene).sum(axis=1)
        beta_cap[:, k] = (D_struct - bonus - log_pi_ratio[k]) / n_nbr

    # When beta_cap < 0 the data and baseline prior already prefer class k
    # over Zero, so no cap is needed and we use the full configured beta.
    effective_real = np.where(beta_cap < 0, beta_cfg, np.minimum(beta_cfg, beta_cap))

    out = np.full((obj.nC, obj.nK), beta_cfg, dtype=np.float32)
    out[:, :-1] = effective_real
    return out


def compute_effective_beta_numba(obj) -> np.ndarray:
    """
    Numba-accelerated equivalent of ``compute_effective_beta``. About 4x
    faster on a 24k-cell run (~1.1s -> ~0.28s on 8 cores). Pre-casts inputs
    to float32 to match the kernel signature, dispatches to the JIT'd
    kernel, then applies the same np.where cap as the numpy version.
    """
    r = np.float32(obj.config['rSpot'])
    spotReg = np.float32(obj.config['SpotReg'])
    n_nbr = np.float32(obj.config['nNeighbors'])
    beta_cfg = obj.config['mrf_beta']

    mu = obj.single_cell.mean_expression_adj.values.astype(np.float32)       # (nG, nK)
    mu_real = mu[:, :-1]                                                     # (nG, nK-1)
    eta_bar = np.asarray(obj.genes.eta_bar, dtype=np.float32)                # (nG,)
    theta_full = np.asarray(obj.cells.theta_bar, dtype=np.float32)           # (nC, nK)
    theta_bar = theta_full[:, :-1]                                          # (nC, nK-1)  theta_{c,k}
    theta_zero = theta_full[:, -1]                                          # (nC,)       theta_{c,Zero}
    geneCount = np.asarray(obj.cells.geneCount, dtype=np.float32)            # (nC, nG)

    log_prior = np.asarray(obj.cellTypes.log_prior, dtype=np.float32)        # (nK,)
    log_pi_ratio = log_prior[:-1] - log_prior[-1]                            # (nK-1,)

    beta_cap = _calc_beta_cap_numba_kernel(mu_real, eta_bar, theta_bar, theta_zero, geneCount,
                                           log_pi_ratio, r, spotReg, n_nbr)

    effective_real = np.where(beta_cap < 0, beta_cfg, np.minimum(beta_cfg, beta_cap))

    out = np.full((obj.nC, obj.nK), beta_cfg, dtype=np.float32)
    out[:, :-1] = effective_real
    return out