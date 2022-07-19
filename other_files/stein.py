import pandas as pd
import numpy as np

"""
implements the Stein covariance estimator
Reference:  https://link.springer.com/content/pdf/10.1007/s00180-016-0672-4.pdf
"""
cov = []
cov.append([[78, 0], [0, 78]])
cov.append([[65, -22], [-22, 80]])
cov = np.array(cov)


def stein(x: np.array, n: np.float) -> np.array:
    """
    x: The covariance matrix
    n: sample size
    """
    p, _p = x.shape
    assert p == _p, "Input is not a squared matrix"
    assert p == 2, "Input must be a 2-by-2 squared matrix"

    [eigval, eigvec] = np.linalg.eig(cov[1])
    # D = np.eye(p)
    # np.fill_diagonal(D, eigval)
    # sigma = eigvec @ D @ eigvec.transpose()

    # p = 2
    # n = 15
    # stein = np.zeros([p, p])

    # Equation 8 from https://link.springer.com/content/pdf/10.1007/s00180-016-0672-4.pdf
    a_0 = n - p + 1 + 2 * eigval[0] + 1 / (eigval[0] - eigval[1])
    a_1 = n - p + 1 + 2 * eigval[1] + 1 / (eigval[1] - eigval[0])

    # Equation 9 from https://link.springer.com/content/pdf/10.1007/s00180-016-0672-4.pdf
    phi = np.zeros([p,p])
    phi[0, 0] = eigvec[0] / a_0
    phi[1, 1] = eigvec[1] / a_1

    sigma = eigvec @ phi @ eigvec.transpose()
    return(sigma)


