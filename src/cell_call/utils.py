import numpy as np
import pandas as pd
import json
import numexpr as ne
import numba as nb
import scipy
import scipy.io as spio
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix, hstack
import os
import config
import urllib, base64
from itertools import accumulate
import gc
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()


def read_image_objects(config):
    img_obj = pd.read_csv(config['expanded_cells'])

    meanCellRadius = np.mean(np.sqrt(img_obj.area / np.pi)) * 0.5
    relCellRadius = np.sqrt(img_obj.area / np.pi) / meanCellRadius

    # append 1 for the misreads
    relCellRadius = np.append(relCellRadius, 1)

    nom = np.exp(-relCellRadius ** 2 / 2) * (1 - np.exp(config['InsideCellBonus'])) + np.exp(config['InsideCellBonus'])
    denom = np.exp(-0.5) * (1 - np.exp(config['InsideCellBonus'])) + np.exp(config['InsideCellBonus'])
    CellAreaFactor = nom / denom

    out = {}
    out['area_factor'] = CellAreaFactor
    out['rel_radius'] = relCellRadius
    out['area'] = np.append(img_obj.area, np.nan)
    out['x'] = np.append(img_obj.x.values, np.nan)
    out['y'] = np.append(img_obj.y.values, np.nan)
    out['cell_id'] = np.append(img_obj.cell_id.values, img_obj.cell_id.shape[0]+1)
    # Last cell is a dummy cell, a super neighbour (ie always a neighbour to any given cell)
    # and will be used to get all the misreads

    return out


def gammaExpectation(rho, beta):
    '''
    :param r:
    :param b:
    :return: Expectetation of a rv X following a Gamma(r,b) distribution with pdf
    f(x;\alpha ,\beta )= \frac{\beta^r}{\Gamma(r)} x^{r-1}e^{-\beta x}
    '''

    # sanity check
    # assert (np.all(rho.coords['cell_id'].data == beta.coords['cell_id'])), 'rho and beta are not aligned'
    # assert (np.all(rho.coords['gene_name'].data == beta.coords['gene_name'])), 'rho and beta are not aligned'
    r = rho[:, :, None]
    b = beta
    gamma = np.empty(b.shape)
    ne.evaluate('r/b', out=gamma)

    # del gamma
    del r
    del b
    gc.collect()
    del gc.garbage[:]
    return gamma


def logGammaExpectation(rho, beta):
    r = rho[:, :, None]
    b = beta
    logb = np.empty(b.shape)
    ne.evaluate("log(b)", out=logb)
    return scipy.special.psi(r) - logb


def _log_gamma(x, y, rho, beta):
    def inner(rho, beta):
        r = rho.data[x, y, None]
        b = beta.data[x, y, :]
        logb = np.empty(b.shape)
        ne.evaluate("log(b)", out=logb)
        log_gamma = scipy.special.psi(r) - logb
        return log_gamma
    return inner



def negBinLoglik(x, r, p):
    '''
    Negative Binomial loglikehood
    :param x:
    :param r:
    :param p:
    :return:
    '''

    # sanity check
    # assert (np.all(da_x.coords['cell_id'].data == da_p.coords['cell_id'])), 'gene counts and beta probabilities are not aligned'
    # assert (np.all(da_x.coords['gene_name'].data == da_p.coords['gene_name'])), 'gene counts and beta probabilities are not aligned'

    contr = np.zeros(p.shape)
    x = x[:, :, None]
    ne.evaluate("x * log(p) + r * log(1 - p)", out=contr)
    return contr


@nb.njit(parallel=True, fastmath=True)
def nb_negBinLoglik(x, r, p):
    '''
    Negative Binomial loglikehood
    :param x:
    :param r:
    :param p:
    :return:
    '''
    out = np.empty(p.shape,p.dtype)

    for i in nb.prange(p.shape[0]):
        for j in range(p.shape[1]):
            if x[i, j, 0] != 0.:
                x_ = x[i, j, 0]
                for k in range(p.shape[2]):
                    out[i, j, k] = x_ * np.log(p[i, j, k]) + r * np.log(1.-p[i, j, k])
            else:
                for k in range(p.shape[2]):
                    out[i, j, k] = r * np.log(1.-p[i, j, k])

    return out



def softmax(X, theta = 1.0, axis = None):
    """
    From https://nolanbconaway.github.io/blog/2017/softmax-numpy
    Compute the softmax of each element along an axis of X.

    Parameters
    ----------
    X: ND-Array. Probably should be floats.
    theta (optional): float parameter, used as a multiplier
        prior to exponentiation. Default = 1.0
    axis (optional): axis to compute values along. Default is the
        first non-singleton axis.

    Returns an array the same size as X. The result will sum to 1
    along the specified axis.
    """

    # make X at least 2d
    y = np.atleast_2d(X)

    # find axis
    if axis is None:
        axis = next(j[0] for j in enumerate(y.shape) if j[1] > 1)

    # multiply y against the theta parameter,
    y = y * float(theta)

    # subtract the max for numerical stability
    y = y - np.expand_dims(np.max(y, axis = axis), axis)

    # exponentiate y
    y = np.exp(y)

    # take the sum along the specified axis
    ax_sum = np.expand_dims(np.sum(y, axis = axis), axis)

    # finally: divide elementwise
    p = y / ax_sum

    # flatten if X was 1D
    if len(X.shape) == 1: p = p.flatten()

    return p


def isConverged(spots, p0, tol):
    p1 = spots.adj_cell_prob
    if p0 is None:
        p0 = np.zeros(p1.shape)
    delta = np.max(np.abs(p1 - p0))
    converged = (delta < tol)
    return converged, delta



