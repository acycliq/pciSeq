import numpy as np
import pandas as pd
import numpy_groupies as npg
from collections.abc import Iterable
from ...src.core.utils.geometry import gaussian_ellipsoid_props, gaussian_contour
import logging

summary_logger = logging.getLogger(__name__)


def cells_summary(cells, spots, genes, is3D):
    '''
    returns a dataframe summarising the main features of each cell, ie gene counts and cell types
    :param spots:
    :return:
    '''
    iCounts = np.argsort(-1 * cells.geneCount, axis=1)
    gene_names = genes.gene_panel[iCounts]
    gene_count = np.take_along_axis(cells.geneCount, iCounts, axis=1)

    iProb = np.argsort(-1 * cells.classProb, axis=1)
    class_names = cells.class_names[iProb]
    class_prob = np.take_along_axis(cells.classProb, iProb, axis=1)

    tol = 0.001

    summary_logger.info('Start collecting data ...')

    isCount_nonZero = [d > tol for d in gene_count]
    name_list = [list(gene_names[i][d]) for (i, d) in enumerate(isCount_nonZero)]
    count_list = [((gene_count[i][d] * 1000).astype(np.int32) / 1000).tolist() for i, d in enumerate(isCount_nonZero)]

    # get spot IDs grouped by cell and gene
    spot_ids = get_contributing_spots_2(cells, spots, genes)

    # reorder each row so spots for genes with higher read counts come first
    spot_ids = np.take_along_axis(spot_ids, iCounts, axis=1)

    # keep only spots with counts above the threshold, and convert each to a plain list
    spot_id_list = [
        [d.tolist() for d in row[mask]]
        for row, mask in zip(spot_ids, isCount_nonZero)
    ]


    isProb_nonZero = [d > tol for d in class_prob]
    class_name_list = [list(class_names[i][d]) for (i, d) in enumerate(isProb_nonZero)]
    prob_list = [((class_prob[i][d] * 1000).astype(np.int32) / 1000).tolist() for i, d in enumerate(isProb_nonZero)]

    contour = []
    for i in range(cells.nC):
        # ea = cells.ellipsoid_attributes[i]
        mu = cells.centroid.iloc[i].tolist()
        cov = cells.cov[i]
        ellipsis = gaussian_contour(mu[:2], cov[:2, :2], 3).astype(np.int64)
        contour.append(ellipsis.tolist())

    df = pd.DataFrame({'Cell_Num': cells.centroid.index.tolist(),
                       'X': ((cells.centroid['x'] * 1000).astype(np.int32) / 1000).tolist(),
                       'Y': ((cells.centroid['y'] * 1000).astype(np.int32) / 1000).tolist(),
                       'Genenames': name_list,
                       'CellGeneCount': count_list,
                       'spot_id': spot_id_list, # the spot ids that when summed-up will generate CellGeneCount
                       'ClassName': class_name_list,
                       'Prob': prob_list,
                       'gaussian_contour': contour
                       })
    if is3D:
        df['sphere_scale'], df['sphere_rotation'] = sphere_props(cells)
        df['Z'] = ((cells.centroid['z'] * 1000).astype(np.int32) / 1000).tolist()
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
                        'spot_id': spots.data.index.tolist(),
                        'x': np.round(spots.data.x.astype('float64'), 3).tolist(),
                        'y': np.round(spots.data.y.astype('float64'), 3).tolist(),
                        'plane_id': spots.data.plane_id.tolist(),
                        'neighbour': max_nbrs.tolist(),
                        'neighbour_array': nbrs.tolist(),
                        'neighbour_prob': p.tolist(),
                        # 'omp_score': ((spots.data.score * 1000).astype(np.int32)/1000).tolist()
                        })
    if is3D:
        out['z'] = np.round(spots.data.z.astype('float64'), 3).tolist()
        out['omp_score'] = np.round(spots.data.score.astype('float64'), 3).tolist()
        out['omp_intensity'] = np.round(spots.data.intensity.astype('float64'), 3).tolist()
        # move column z after x, y
        z_pos = out.columns.get_loc('y') + 1
        out.insert(z_pos, 'z', out.pop('z'))

    return out


def collect_data(cells, spots, genes, is3D):
    '''
    Collects data for the viewer
    :param cells:
    :param spots:
    :return:
    '''
    cell_df = cells_summary(cells, spots, genes, is3D)
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


def get_contributing_spots(spot_ids, gene_id, parent_cell_id, parent_cell_prob, tol=0.001):
    """
    Get only the spot lists for each (cell, gene) combination.
    Stripped down version of aggregate_cell_gene_with_spots.

    Args:
        spot_ids: array of shape (n_spots,) - spot identifiers
        gene_id: array of shape (n_spots,) - gene category for each spot
        parent_cell_id: array of shape (n_spots, n_assignments) - cell labels
        parent_cell_prob: array of shape (n_spots, n_assignments) - probabilities

    Returns:
        spot_lists: array of shape (n_cells, n_genes) - lists of contributing spot_ids (prob > 0.001)
        cell_ids: array of unique cell IDs (row labels)
        gene_ids: array of unique gene IDs (column labels)
    """

    # Get unique cell and gene IDs
    unique_cells = np.unique(parent_cell_id.ravel())
    unique_genes = np.unique(gene_id)

    # Create mapping dictionaries
    cell_to_idx = {cell: idx for idx, cell in enumerate(unique_cells)}
    gene_to_idx = {gene: idx for idx, gene in enumerate(unique_genes)}

    n_cells = len(unique_cells)
    n_genes = len(unique_genes)

    # Expand arrays to match flattened structure
    spot_expanded = np.repeat(spot_ids, parent_cell_id.shape[1])
    gene_expanded = np.repeat(gene_id, parent_cell_id.shape[1])

    # Flatten the 2D arrays
    cells_flat = parent_cell_id.ravel()
    probs_flat = parent_cell_prob.ravel()

    # Filter valid entries (probabilities > 0.001 and valid cells)
    valid_mask = (probs_flat > tol) & np.isin(cells_flat, unique_cells)

    spot_valid = spot_expanded[valid_mask]
    gene_valid = gene_expanded[valid_mask]
    cells_valid = cells_flat[valid_mask]

    # Convert to indices
    gene_indices = np.array([gene_to_idx[g] for g in gene_valid])
    cell_indices = np.array([cell_to_idx[c] for c in cells_valid])

    # Create linear group indices
    group_indices = cell_indices * n_genes + gene_indices

    # Initialize spot lists arrays - REMOVED probability aggregation
    spot_lists_flat = np.empty(n_cells * n_genes, dtype=object)
    for i in range(n_cells * n_genes):
        spot_lists_flat[i] = []

    # Collect spots for each group - ONLY spot collection, no probability aggregation
    for i, group_id in enumerate(group_indices):
        spot_lists_flat[group_id].append(spot_valid[i].item())

    # Reshape results - REMOVED cell_gene_matrix and prob_lists
    spot_lists = spot_lists_flat.reshape(n_cells, n_genes)

    return spot_lists


def get_contributing_spots_2(cells, spots, genes):
    nN = spots.parent_cell_id.shape[1]
    cell_ids = spots.parent_cell_id.ravel()
    gene_ids = np.tile(spots.gene_id, (nN, 1)).T.ravel()
    group_idx =  np.vstack((cell_ids, gene_ids))

    spot_ids = np.tile(spots.data.index.values, (nN, 1)).T
    spot_ids = spot_ids.ravel()

    agg = npg.aggregate_np(group_idx, spot_ids, size=(cells.nC, genes.nG), func=list, fill_value=[], dtype=object)

    # # assert cells.geneCount.shape == agg.shape
    # out = []
    # for i, counts in enumerate(agg):
    #     counts = [list(d) for d in counts]
    #     out.append(counts)
    # # #     spot_id_list = agg[i][mask].tolist()
    # #     spot_id_list = agg[i].tolist()
    # #     out.append(counts)
    return agg