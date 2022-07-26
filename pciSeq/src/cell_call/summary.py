import numpy as np
import pandas as pd
from pciSeq.src.cell_call.utils import gaussian_ellipsoid, gaussian_ellipsoid_props
from pciSeq.src.cell_call.log_config import logger


def _iss_summary(cells, genes, single_cell, cfg):
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
    cov_list = []
    if cfg['is_3D']:
        sphere_scale = []
        sphere_rotation = []
        cov_list = []
        for i in range(cells.nC):
            cov = cells.cov[i]
            _sphere_scale, _sphere_rotation = gaussian_ellipsoid_props(cov, 3)
            sphere_scale.append(_sphere_scale.tolist())
            sphere_rotation.append(_sphere_rotation)
            cov_list.append(cov.flatten().tolist())

    if not cfg['is_3D'] and cfg['relax_segmentation']:
        for i in range(cells.nC):
            # ea = cells.ellipsoid_attributes[i]
            mu = cells.centroid.iloc[i].tolist()
            cov = cells.cov[i]
            ellipsis = gaussian_ellipsoid(mu[:2], cov[:2, :2], 3).astype(np.int)
            ellipsoid_border.append(ellipsis.tolist())

    iss_df = pd.DataFrame({'Cell_Num': cells.cell_props['cell_label'].tolist(),
                           'X': cells.cell_props['x'].tolist(),
                           'Y': cells.cell_props['y'].tolist(),
                           'Genenames': name_list,
                           'CellGeneCount': count_list,
                           'ClassName': class_name_list,
                           'Prob': prob_list
                            },
                           columns=['Cell_Num', 'X', 'Y', 'Genenames', 'CellGeneCount', 'ClassName', 'Prob']
                           )
    if np.any(cells.cell_props['z']):
        iss_df = iss_df.assign(Z=cells.cell_props['z'].tolist())

    if len(sphere_scale) > 0:
        iss_df = iss_df.assign(sphere_scale = sphere_scale)

    if len(sphere_rotation) > 0:
        iss_df = iss_df.assign(sphere_rotation = sphere_rotation)

    if len(cov_list) > 0:
        iss_df = iss_df.assign(cov = cov_list)

    if len(ellipsoid_border) > 0:
        iss_df = iss_df.assign(ellipsoid_border = ellipsoid_border)

    iss_df.set_index(['Cell_Num'])

    # Ignore the first row. It is the pseudocell to keep the misreads (ie the background)
    iss_df = iss_df[1:]
    logger.info(' Data collected!')

    return iss_df


def _summary(spots):
    # check for duplicates (ie spots with the same coordinates with or without the same gene name).
    # is_duplicate = spots.data.duplicated(subset=['x', 'y'])

    num_rows = spots.data.shape[0]

    cell_prob = spots.parent_cell_prob
    neighbors = spots.parent_cell_id
    p = [cell_prob[i, :].tolist() for i in range(num_rows)]
    nbrs = [neighbors[i, :].tolist() for i in range(num_rows)]
    max_nbrs = [neighbors[i, idx].tolist() for i in range(num_rows) for idx in [np.argmax(cell_prob[i, :])]]

    out = pd.DataFrame({'Gene': spots.data.Gene.tolist(),
                        'Gene_id': spots.gene_id.tolist(),
                        'x': spots.data.x.tolist(),
                        'y': spots.data.y.tolist(),
                        'neighbour': max_nbrs,
                        'neighbour_array': nbrs,
                        'neighbour_prob': p})

    if np.any(spots.data.z):
        out = out.assign(z=spots.data.z.tolist())

    return out


def collect_data(cells, spots, genes, single_cell, cfg):
    '''
    Collects data for the viewer
    :param cells:
    :param spots:
    :return:
    '''
    iss_df = _iss_summary(cells, genes, single_cell, cfg)
    gene_df = _summary(spots)
    return iss_df, gene_df
