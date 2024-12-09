import numpy as np
import pandas as pd
from ...src.core.utils.geometry import gaussian_ellipsoid_props, gaussian_contour
import logging

summary_logger = logging.getLogger(__name__)


def cells_summary(cells, genes, single_cell, is3D):
    '''
    returns a dataframe summarising the main features of each cell, ie gene counts and cell types
    :param spots:
    :return:
    '''
    iCounts = np.argsort(-1 * cells.geneCount, axis=1)
    gene_names = genes.gene_panel[iCounts]
    gene_count = np.take_along_axis(cells.geneCount, iCounts, axis=1)

    iProb = np.argsort(-1 * cells.classProb, axis=1)
    class_names = single_cell.classes[iProb]
    class_prob = np.take_along_axis(cells.classProb, iProb, axis=1)

    tol = 0.001

    summary_logger.info('Start collecting data ...')

    isCount_nonZero = [d > tol for d in gene_count]
    name_list = [list(gene_names[i][d]) for (i, d) in enumerate(isCount_nonZero)]
    count_list = [list(gene_count[i][d].round(3)) for (i, d) in enumerate(isCount_nonZero)]

    isProb_nonZero = [d > tol for d in class_prob]
    class_name_list = [list(class_names[i][d]) for (i, d) in enumerate(isProb_nonZero)]
    prob_list = [list(class_prob[i][d].round(3)) for (i, d) in enumerate(isProb_nonZero)]

    contour = []
    for i in range(cells.nC):
        # ea = cells.ellipsoid_attributes[i]
        mu = cells.centroid.iloc[i].tolist()
        cov = cells.cov[i]
        ellipsis = gaussian_contour(mu[:2], cov[:2, :2], 3).astype(np.int64)
        contour.append(ellipsis.tolist())

    df = pd.DataFrame({'Cell_Num': cells.centroid.index.tolist(),
                       'X': cells.centroid['x'].round(3).tolist(),
                       'Y': cells.centroid['y'].round(3).tolist(),
                       'Genenames': name_list,
                       'CellGeneCount': count_list,
                       'ClassName': class_name_list,
                       'Prob': prob_list,
                       'gaussian_contour': contour
                       })
    if is3D:
        df['sphere_scale'], df['sphere_rotation'] = sphere_props(cells)
        df['Z'] = cells.centroid['z'].tolist()
        # move column Z after X, Y
        df.insert(3, 'Z', df.pop('Z'))

    df.set_index(['Cell_Num'])

    # Ignore the first row. It is the pseudocell to keep the misreads (ie the background)
    df = df[1:]
    summary_logger.info('Data collected!')
    return df


def spots_summary(spots, is3D):
    # check for duplicates (ie spots with the same coordinates with or without the same gene name).
    # is_duplicate = spots.data.duplicated(subset=['x', 'y'])

    idx = np.argsort(-1 * spots.parent_cell_prob, axis=1)
    p = np.take_along_axis(spots.parent_cell_prob, idx, axis=1).round(3)
    nbrs = np.take_along_axis(spots.parent_cell_id, idx, axis=1)
    max_nbrs = nbrs[:, 0]

    out = pd.DataFrame({'gene_name': spots.data.gene_name.tolist(),
                        'gene_id': spots.gene_id.tolist(),
                        'x': spots.data.x.tolist(),
                        'y': spots.data.y.tolist(),
                        'neighbour': max_nbrs.tolist(),
                        'neighbour_array': nbrs.tolist(),
                        'neighbour_prob': p.tolist()
                        })
    if is3D:
        out['z'] = spots.data.z.tolist()
        # move column z after x, y
        out.insert(4, 'z', out.pop('z'))

    return out


def collect_data(cells, spots, genes, single_cell, is3D):
    '''
    Collects data for the viewer
    :param cells:
    :param spots:
    :return:
    '''
    cell_df = cells_summary(cells, genes, single_cell, is3D)
    gene_df = spots_summary(spots, is3D)
    return cell_df, gene_df


def sphere_props(cells):
    sphere_scale = []
    sphere_rotation = []
    for i in range(cells.nC):
        cov = cells.cov[i]
        scale, rotation = gaussian_ellipsoid_props(cov)
        sphere_scale.append(scale)
        sphere_rotation.append(rotation)
    return sphere_scale, sphere_rotation
