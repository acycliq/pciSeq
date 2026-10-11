import numpy as np
import pandas as pd
import numpy_groupies as npg
from collections.abc import Iterable
from ...src.core.utils.geometry import gaussian_ellipsoid_props, gaussian_contour
import logging

logger = logging.getLogger(__name__)


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

    logger.info('Start collecting data ...')

    isCount_nonZero = [d > tol for d in gene_count]
    name_list = [list(gene_names[i][d]) for (i, d) in enumerate(isCount_nonZero)]
    count_list = [((gene_count[i][d] * 1000).astype(np.int32) / 1000).tolist() for i, d in enumerate(isCount_nonZero)]

    # get spot IDs grouped by cell and gene, keeping a spot only if its probability for
    # that cell is above spot_tol.
    #
    # tol / 10 and not tol: tol is the cut-off for CellGeneCount, and a count just above
    # it can be the sum of several spots that are each below it. Cutting the spots at tol
    # too would show such a gene with a count and no spots. Not a guarantee, but any
    # mismatch left can only be around tol.
    spot_tol = tol / 10
    spot_ids = get_contributing_spots_2(cells, spots, genes, spot_tol)

    # reorder each row so spots for genes with higher read counts come first
    spot_ids = np.take_along_axis(spot_ids, iCounts, axis=1)

    # keep only the genes with counts above the threshold, and convert each to a plain
    # list. np.asarray first: a gene with no spot above spot_tol comes back as the empty
    # python list that aggregate fills in, and that has no tolist of its own.
    spot_id_list = [
        [np.asarray(d).tolist() for d in row[mask]]
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
    # the segmentation centroid as well, X0, Y0, Z0. X, Y, Z is the fitted centre
    # and the two drift apart when cell_centroid_prior is small, so keep where
    # the cell started
    m0 = cells.ini_centroids()
    if is3D:
        df['sphere_scale'], df['sphere_rotation'] = sphere_props(cells)
        df['Z'] = ((cells.centroid['z'] * 1000).astype(np.int32) / 1000).tolist()
        df['Z0'] = ((m0['z'] * 1000).astype(np.int32) / 1000).tolist()
    else:
        # 2D has no z to report, but the viewer's arrow export asks for the
        # column whatever the data is, so give it zeros. Same idea as z_plane
        # on the spots side. Without this cellData_to_arrow refuses to write
        # and a 2D run comes out with no cells for the viewer.
        df['Z'] = np.float32(0.0)
        df['Z0'] = np.float32(0.0)
    # move column Z after X, Y, then X0, Y0, Z0 after that
    df.insert(3, 'Z', df.pop('Z'))
    df.insert(4, 'X0', ((m0['x'] * 1000).astype(np.int32) / 1000).tolist())
    df.insert(5, 'Y0', ((m0['y'] * 1000).astype(np.int32) / 1000).tolist())
    df.insert(6, 'Z0', df.pop('Z0'))

    df.set_index(['Cell_Num'])

    # Ignore the first row. It is the pseudocell to keep the misreads (ie the background)
    df = df[1:]
    logger.info('Data collected!')
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
    # These three used to be written only for 3D, but the viewer reads all of
    # them whatever the data is, so a 2D run left it with nothing to draw. In 2D
    # z is zero for every spot anyway, and score and intensity are filled in by
    # the validator when the caller has none.
    out['z'] = np.round(spots.data.z.astype('float64'), 3).tolist()
    out['omp_score'] = np.round(spots.data.score.astype('float64'), 3).tolist()
    out['omp_intensity'] = np.round(spots.data.intensity.astype('float64'), 3).tolist()
    # which cell the spot physically sits in, straight off the label image, 0 if none.
    # Not the same as 'neighbour': that is the cell the model decided the spot came
    # from, this is just where the segmentation put it. The model only ever looks at
    # it through InsideCellBonus, which is normally 0. Named after inside_cell in
    # label_processing, which is what works it out.
    out['inside_cell'] = spots.data.label.tolist()
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


def get_contributing_spots_2(cells, spots, genes, tol=0.0):
    """For every cell and gene, the ids of the spots of that gene that belong to the cell.

    A spot counts as belonging to a cell when its probability for that cell is above
    tol. With tol at 0 every candidate spot is listed, which is what this did before.
    Returns an object array, cells by genes, with a list of spot ids in each entry.
    """
    nN = spots.parent_cell_id.shape[1]
    cell_ids = spots.parent_cell_id.ravel()
    gene_ids = np.tile(spots.gene_id, (nN, 1)).T.ravel()

    spot_ids = np.tile(spots.data.index.values, (nN, 1)).T
    spot_ids = spot_ids.ravel()

    # parent_cell_prob lines up with parent_cell_id, one column per candidate cell
    keep = spots.parent_cell_prob.ravel() > tol
    group_idx = np.vstack((cell_ids[keep], gene_ids[keep]))

    agg = npg.aggregate_np(group_idx, spot_ids[keep], size=(cells.nC, genes.nG), func=list, fill_value=[], dtype=object)

    # # assert cells.geneCount.shape == agg.shape
    # out = []
    # for i, counts in enumerate(agg):
    #     counts = [list(d) for d in counts]
    #     out.append(counts)
    # # #     spot_id_list = agg[i][mask].tolist()
    # #     spot_id_list = agg[i].tolist()
    # #     out.append(counts)
    return agg