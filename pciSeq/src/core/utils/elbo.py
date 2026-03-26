import logging

import numpy as np

logger = logging.getLogger(__name__)


def calc_elbo(obj):
    E_log_p = expected_log_joint(obj)
    entropy = total_entropy(obj)
    return E_log_p + entropy


def total_entropy(obj):
    """
    Computes the total entropy of the variational distribution: H[q].
    = H[q(z)] + H[q(zeta)] + H[q(gamma)] + H[q(eta)]
    No entropy for theta (point estimate).
    """
    H_z = categorical_entropy(obj.spots.parent_cell_prob)
    H_zeta = categorical_entropy(obj.cells.classProb)
    H_gamma = entropy_gamma(obj)
    H_eta = entropy_eta(obj)
    H_pi = entropy_pi(obj)

    nS = obj.spots.parent_cell_prob.shape[0]
    nN = obj.spots.parent_cell_prob.shape[1]
    nC = obj.cells.classProb.shape[0]
    nK = obj.cells.classProb.shape[1]
    logger.info('Entropy: H[z]=%.2f (%.4f/spot, max=%.2f) | H[zeta]=%.2f (%.4f/cell, max=%.2f) | H[gamma]=%.2f | H[eta]=%.2f | H[pi]=%.2f',
                H_z, H_z / nS, np.log(nN),
                H_zeta, H_zeta / nC, np.log(nK),
                H_gamma, H_eta, H_pi)

    return H_z + H_zeta + H_gamma + H_eta + H_pi


def expected_log_joint(obj):
    """
    Computes E_q[log p(data, latent variables)] = likelihood + all priors.
    """
    log_lik = poisson_process_loglikelihood(obj)
    log_prior_zeta = zeta_prior(obj)
    log_prior_mrf = mrf_prior(obj)
    log_prior_theta = theta_prior(obj)
    log_prior_gamma = gamma_prior(obj)
    log_prior_eta = eta_prior(obj)
    log_prior_pi = pi_prior(obj)

    return log_lik + log_prior_zeta + log_prior_mrf + log_prior_theta + log_prior_gamma + log_prior_eta + log_prior_pi

def poisson_process_loglikelihood(obj):
    """
    Computes E_q[log p(spots | z, zeta, theta, gamma, eta)].

    The Poisson process log-likelihood has two parts:

    term1 (integrated intensity): The expected total rate over all cells, penalizing
        the model for predicting spots that were not observed.
        = -sum_{c,g,k} q(zeta_ck) * theta_ck * mu_gk * A_c * E[gamma_cgk] * E[eta_g]

    term2 (observed spots): For each observed spot, the log-rate at its location,
        marginalized over spot-to-cell assignments q(z) and cell types q(zeta).
        Split into:
        - term2_bg: spots assigned to background, weighted by q(z_s = background)
        - term2_sig: spots assigned to cells, summing over neighbor cells and classes
    """
    zeta_ck = obj.cells.classProb
    theta_ck = obj.cells.theta_bar
    mu_gk = obj.single_cell.mean_expression_adj.values + obj.config['SpotReg']
    A_c = obj.cells.ini_cell_props['area_factor']
    gamma_cgk = obj.spots.gamma_bar.compute()
    eta_g = obj.genes.eta_bar

    # term1: -integral of expected intensity over all cells
    term1 = -np.einsum('ck, ck, gk, c, cgk, g -> ', zeta_ck, theta_ck, mu_gk, A_c, gamma_cgk, eta_g)

    # term2: sum over observed spots of log(rate at spot location)
    spots = obj.spots
    cfg = obj.config
    nN = cfg['nNeighbors'] + 1
    g_s = spots.gene_id                    # (nS,) gene index per spot
    q_z = spots.parent_cell_prob           # (nS, nN) — q(z_s)
    cell_ids = spots.parent_cell_id        # (nS, nN)
    log_mu = np.log(mu_gk)                 # (nG, nK)
    log_theta = np.log(theta_ck)           # (nC, nK) — point estimate, so E[log] = log
    E_log_gamma = spots.log_gamma_bar      # (nC, nG, nK) — E_q[log gamma] = psi(a) - log(b)
    if hasattr(E_log_gamma, 'compute'):
        E_log_gamma = E_log_gamma.compute()
    E_log_eta = obj.genes.logeta_bar       # (nG,) — E_q[log eta] = psi(a) - log(b)
    mvn_loglik = spots.mvn_loglik_arr      # (nS, nN) — spatial: multivariate normal log-pdf
    bonus = spots.bonus_mask * cfg['InsideCellBonus']

    # Background: q(z_s = bg) * E[log rho_g]
    log_rho = obj.genes.log_rho_bar[g_s]
    term2_bg = np.sum(q_z[:, -1] * log_rho)

    # Signal: for each neighbor cell, sum over classes weighted by q(zeta)
    log_mu_s = log_mu[g_s, :]              # (nS, nK) — precompute, same for all neighbors
    term2_sig = 0.0
    for n in range(nN - 1):
        c_n = cell_ids[:, n]

        # sum_k q(zeta_{c_n, k}) * [log mu_{g_s,k} + log theta_{c_n,k} + E[log gamma_{c_n,g_s,k}]]
        class_term = np.sum(
            zeta_ck[c_n, :] * (log_theta[c_n, :] + log_mu_s + E_log_gamma[c_n, g_s, :]),
            axis=1
        )  # (nS,)

        # weight by q(z_s = c_n) and add spatial + gene efficiency terms
        term2_sig += np.sum(q_z[:, n] * (mvn_loglik[:, n] + bonus[:, n] + class_term + E_log_eta[g_s]))

    return term1 + term2_bg + term2_sig


def zeta_prior(obj):
    """
    Computes E_q[log p(zeta | pi)].

    Each cell's type is drawn from a Categorical with prior probabilities pi_k.
    = sum_{c,k} q(zeta_ck) * log(pi_k)
    """
    zeta_ck = obj.cells.classProb       # (nC, nK)
    log_pi = obj.cellTypes.log_prior    # (nK,)
    return np.sum(zeta_ck * log_pi)


def mrf_prior(obj):
    """
    Computes the MRF (Potts model) contribution to the expected log-joint.

    Encourages neighboring cells to have the same type.
    = beta * sum_{c,k} q(zeta_ck) * sum_{c' in neighbors(c)} q(zeta_{c',k})
    """
    zeta_ck = obj.cells.classProb                           # (nC, nK)
    beta = obj.config['mrf_beta']
    mrf_val = zeta_ck[obj.cells.nbrs].sum(axis=1)           # (nC, nK)
    return beta * np.sum(zeta_ck * mrf_val)


def theta_prior(obj):
    """
    Computes E_q[log p(theta | zeta)].

    Theta has a Gamma(r, r) prior, conditional on zeta, so weighted by q(zeta_ck).
    = sum_{c,k} q(zeta_ck) * [r*log(r) - gammaln(r) + (r-1)*log(theta_ck) - r*theta_ck]

    Theta is a point estimate, so E[log theta] = log(theta_bar).
    """
    from scipy.special import gammaln

    zeta_ck = obj.cells.classProb       # (nC, nK)
    theta_ck = obj.cells.theta_bar      # (nC, nK)
    r = obj.config['rTheta']

    # log p(theta) = r*log(r) - gammaln(r) + (r-1)*log(theta) - r*theta
    log_norm = r * np.log(r) - gammaln(r)                       # normalising constant of the Gamma pdf
    log_pdf = log_norm + (r - 1) * np.log(theta_ck) - r * theta_ck  # log Gamma pdf evaluated at theta_ck

    return np.sum(zeta_ck * log_pdf)


def gamma_prior(obj):
    """
    Computes E_q[log p(gamma | zeta)].

    Gamma has a Gamma(r, r) prior, conditional on zeta, so weighted by q(zeta_ck).
    = sum_{c,g,k} q(zeta_ck) * [r*log(r) - gammaln(r) + (r-1)*E[log gamma_cgk] - r*E[gamma_cgk]]

    Gamma has a full variational distribution, so we use:
        E[log gamma] = psi(alpha) - log(beta)  (stored as log_gamma_bar)
        E[gamma] = alpha / beta                 (stored as gamma_bar)
    """
    from scipy.special import gammaln

    zeta_ck = obj.cells.classProb               # (nC, nK)
    r = obj.config['rSpot']

    E_gamma = obj.spots.gamma_bar               # (nC, nG, nK)
    if hasattr(E_gamma, 'compute'):
        E_gamma = E_gamma.compute()
    E_log_gamma = obj.spots.log_gamma_bar       # (nC, nG, nK)
    if hasattr(E_log_gamma, 'compute'):
        E_log_gamma = E_log_gamma.compute()

    # log p(gamma) = r*log(r) - gammaln(r) + (r-1)*E[log gamma] - r*E[gamma]
    log_norm = r * np.log(r) - gammaln(r)       # normalising constant of the Gamma pdf
    log_pdf = log_norm + (r - 1) * E_log_gamma - r * E_gamma  # (nC, nG, nK)

    # weight by q(zeta_ck), broadcast over genes: (nC, 1, nK)
    zeta_weight = zeta_ck[:, None, :]
    return np.sum(zeta_weight * log_pdf)


def eta_prior(obj):
    """
    Computes E_q[log p(eta)].

    Eta has a Gamma(r, r) prior, independent of zeta (no weighting needed).
    = sum_g [r*log(r) - gammaln(r) + (r-1)*E[log eta_g] - r*E[eta_g]]
    """
    from scipy.special import gammaln

    r = obj.config['rGene']
    E_eta = obj.genes.eta_bar           # (nG,)
    E_log_eta = obj.genes.logeta_bar    # (nG,)

    # log p(eta) = r*log(r) - gammaln(r) + (r-1)*E[log eta] - r*E[eta]
    log_norm = r * np.log(r) - gammaln(r)                   # normalising constant of the Gamma pdf
    log_pdf = log_norm + (r - 1) * E_log_eta - r * E_eta  # log Gamma pdf evaluated at each gene

    return np.sum(log_pdf)


def pi_prior(obj):
    """
    Computes E_q[log p(pi | alpha_0)] for the real classes only.

    Only active in 'weighted' mode where pi is a Dirichlet random variable.
    In 'uniform' mode pi is fixed, so this term is zero.
    Zero class is excluded (its weight is fixed, not part of the Dirichlet).

    E_q[log p(pi | alpha_0)] = -log B(alpha_0) + sum_k (alpha_0k - 1) * E[log pi_k]
    where E[log pi_k] = psi(alpha_post_k) - psi(sum(alpha_post))
    """
    if obj.config['cell_type_prior'] != 'weighted' and not obj.single_cell.isMissing:
        return 0.0

    from scipy.special import gammaln, psi

    # alpha is (nK-1,) — real classes only, prior alpha_0 = ones
    alpha_0 = np.ones(obj.cellTypes.nK - 1, dtype=np.float32)
    alpha_post = obj.cellTypes.alpha

    # E[log pi_k] under the Dirichlet posterior (real classes only)
    E_log_pi = psi(alpha_post) - psi(alpha_post.sum())

    # -log B(alpha_0) = gammaln(sum(alpha_0)) - sum(gammaln(alpha_0k))
    log_norm = gammaln(alpha_0.sum()) - np.sum(gammaln(alpha_0))

    return log_norm + np.sum((alpha_0 - 1) * E_log_pi)


def entropy_pi(obj):
    """
    Entropy of q(pi) = Dirichlet(alpha_post).

    Only active in 'weighted' mode. In 'uniform' mode pi is fixed (no entropy).

    H[Dirichlet(alpha)] = log B(alpha) + (sum(alpha) - K) * psi(sum(alpha))
                          - sum_k (alpha_k - 1) * psi(alpha_k)
    """
    if obj.config['cell_type_prior'] != 'weighted' and not obj.single_cell.isMissing:
        return 0.0

    from scipy.special import gammaln, psi

    alpha = obj.cellTypes.alpha
    K = len(alpha)
    alpha_sum = alpha.sum()

    # log B(alpha) = sum(gammaln(alpha_k)) - gammaln(sum(alpha))
    log_B = np.sum(gammaln(alpha)) - gammaln(alpha_sum)

    return log_B + (alpha_sum - K) * psi(alpha_sum) - np.sum((alpha - 1) * psi(alpha))


def categorical_entropy(probs):
    """
    Entropy of a categorical distribution: H = -sum(p * log(p)).
    Handles zeros safely (0 * log(0) = 0).

    Args:
        probs: array of probabilities, e.g. (nS, nN) or (nC, nK)
    """
    safe_probs = np.where(probs > 0, probs, 1.0)
    return -np.sum(probs * np.log(safe_probs))


def entropy_gamma(obj):
    """
    Entropy of q(gamma), weighted by q(zeta_ck) since gamma is conditional on zeta.

    H = sum_{c,g,k} q(zeta_ck) * [alpha - log(beta) + gammaln(alpha) + (1-alpha)*psi(alpha)]

    where alpha = _post_shape and beta = _post_rate of the variational Gamma posterior.
    """
    from scipy.special import gammaln, psi

    zeta_ck = obj.cells.classProb           # (nC, nK)
    alpha = obj.spots._post_shape           # (nC, nG) or (nC, nG, nK)
    beta = obj.spots._post_rate             # (nC, nG, nK)

    # alpha may be (nC, nG) — broadcast over K
    if alpha.ndim == 2:
        alpha = alpha[:, :, None]

    # H[Gamma(alpha, beta)] = alpha - log(beta) + gammaln(alpha) + (1 - alpha) * psi(alpha)
    entropy = alpha - np.log(beta) + gammaln(alpha) + (1 - alpha) * psi(alpha)

    # weight by q(zeta_ck), summing over cells, genes and classes
    return np.einsum('ck, cgk ->', zeta_ck, entropy)


def entropy_eta(obj):
    """
    Entropy of q(eta). Not conditional on zeta, so no weighting.

    H = sum_g [alpha - log(beta) + gammaln(alpha) + (1-alpha)*psi(alpha)]
    """
    from scipy.special import gammaln, psi

    alpha = obj.genes._post_shape           # (nG,)
    beta = obj.genes._post_rate             # (nG,)

    # H[Gamma(alpha, beta)] = alpha - log(beta) + gammaln(alpha) + (1 - alpha) * psi(alpha)
    entropy = alpha - np.log(beta) + gammaln(alpha) + (1 - alpha) * psi(alpha)

    return np.sum(entropy)
