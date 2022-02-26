import os
import numpy as np
import pandas as pd
from pciSeq.src.cell_call.utils import gaussian_ellipsoid, gaussian_ellipsoid_props
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger(__name__)


def _iss_summary(cells, genes, single_cell, ppm):
    '''
    returns a dataframe summarising the main features of each cell, ie gene counts and cell types
    :param spots:
    :return:
    '''
    x = cells.cell_props['x']
    y = cells.cell_props['y']
    cell_id = cells.cell_props['cell_label']

    gene_count = cells.geneCount
    class_prob = cells.classProb
    gene_names = genes.gene_panel
    class_names = single_cell.classes

    tol = 0.001

    logger.info(' Start collecting data ...')
    N = len(cell_id)
    isCount_nonZero = [gene_count[n, :] > tol for n in range(N)]
    name_list = [gene_names[isCount_nonZero[n]].tolist() for n in range(N)]
    count_list = [gene_count[n, isCount_nonZero[n]].tolist() for n in range(N)]

    isProb_nonZero = [class_prob[n, :] > tol for n in range(N)]
    class_name_list = [class_names[isProb_nonZero[n]].tolist() for n in range(N)]
    prob_list = [class_prob[n, isProb_nonZero[n]].tolist() for n in range(N)]

    ellipsoid_border = []
    sphere_scale = []
    sphere_rotation = []
    for i in range(cells.nC):
        ea = cells.ellipsoid_attributes[i]
        cov = cells.cov[i] * ppm
        mu = cells.centroid.values[i] * ppm
        ellipsis = gaussian_ellipsoid(mu, cov, 3).astype(np.int)
        ellipsoid_border.append(ellipsis.tolist())
        _sphere_scale, _sphere_rotation = gaussian_ellipsoid_props(cov, 3)
        sphere_scale.append(_sphere_scale.tolist())
        sphere_rotation.append(_sphere_rotation)

    iss_df = pd.DataFrame({'Cell_Num': cells.cell_props['cell_label'].tolist(),
                           'X': (ppm * cells.centroid.x).tolist(),
                           'Y': (ppm * cells.centroid.y).tolist(),
                           'Z': (ppm * cells.centroid.z).tolist(),
                           'X_0': (ppm * cells.cell_props['x']).tolist(),
                           'Y_0': (ppm * cells.cell_props['y']).tolist(),
                           'Genenames': name_list,
                           'CellGeneCount': count_list,
                           'ClassName': class_name_list,
                           'Prob': prob_list,
                           'rho': cells.corr.tolist(),
                           'sigma_x': cells.sigma_x.tolist(),
                           'sigma_y': cells.sigma_y.tolist(),
                           'ellipsoid_border': ellipsoid_border,
                           'sphere_scale': sphere_scale,
                           'sphere_rotation': sphere_rotation,
                            },
                          columns=['Cell_Num', 'X', 'Y', 'Z', 'X_0', 'Y_0', 'Genenames', 'CellGeneCount', 'ClassName',
                                   'Prob', 'rho', 'sigma_x', 'sigma_y', 'ellipsoid_border', 'sphere_position',
                                   'sphere_scale', 'sphere_rotation']
                          )
    iss_df.set_index(['Cell_Num'])

    # Ignore the first row. It is the pseudocell to keep the misreads (ie the background)
    iss_df = iss_df[1:]
    logger.info(' Data collected!')

    return iss_df


def _summary(spots, ppm):
    # check for duplicates (ie spots with the same coordinates with or without the same gene name).
    # is_duplicate = spots.data.duplicated(subset=['x', 'y'])

    num_rows = spots.data.shape[0]

    cell_prob = spots.parent_cell_prob
    neighbors = spots.parent_cell_id
    p = [cell_prob[i, :].tolist() for i in range(num_rows)]
    nbrs = [neighbors[i, :].tolist() for i in range(num_rows)]
    max_nbrs = [neighbors[i, idx].tolist() for i in range(num_rows) for idx in [np.argmax(cell_prob[i, :])]]

    out = pd.DataFrame({'Gene': spots.data.gene_name.tolist(),
                        'Gene_id': spots.gene_id.tolist(),
                        'x': (ppm * spots.data.x).tolist(),
                        'y': (ppm * spots.data.y).tolist(),
                        'z': (ppm * spots.data.z).tolist(),
                        'neighbour': max_nbrs,
                        'neighbour_array': nbrs,
                        'neighbour_prob': p})

    return out


def collect_data(cells, spots, genes, single_cell, ppm):
    '''
    Collects data for the viewer
    :param cells:
    :param spots:
    :return:
    '''
    iss_df = _iss_summary(cells, genes, single_cell, ppm)
    gene_df = _summary(spots, ppm)
    return iss_df, gene_df
