import numpy as np
import pandas as pd
import logging

summary_logger = logging.getLogger(__name__)


def _iss_summary(cells, genes, single_cell):
    '''
    returns a dataframe summarising the main features of each cell, ie gene counts and cell types
    :param spots:
    :return:
    '''
    iCounts = [np.argsort(-1 * d) for d in cells.geneCount]
    iProb = [np.argsort(-1 * d) for d in cells.classProb]

    gene_names = [genes.gene_panel[d] for d in iCounts]
    gene_count = [cells.geneCount[i][d] for i, d in enumerate(iCounts)]

    class_names = [single_cell.classes[d] for d in iProb]
    class_prob = [cells.classProb[i][d] for i, d in enumerate(iProb)]

    tol = 0.001

    summary_logger.info('Start collecting data ...')
    isCount_nonZero = [d > tol for d in gene_count]
    name_list = [list(np.array(gene_names[i])[d]) for (i, d) in enumerate(isCount_nonZero)]
    count_list = [list(np.array(gene_count[i])[d]) for (i, d) in enumerate(isCount_nonZero)]

    isProb_nonZero = [d > tol for d in class_prob]
    class_name_list = [list(np.array(class_names[i])[d]) for (i, d) in enumerate(isProb_nonZero)]
    prob_list = [list(np.array(class_prob[i])[d]) for (i, d) in enumerate(isProb_nonZero)]

    iss_df = pd.DataFrame({'Cell_Num': cells.centroid.index.tolist(),
                           'X': cells.centroid['x'].tolist(),
                           'Y': cells.centroid['y'].tolist(),
                           'Genenames': name_list,
                           'CellGeneCount': count_list,
                           'ClassName': class_name_list,
                           'Prob': prob_list
                            },
                           columns=['Cell_Num', 'X', 'Y', 'Genenames', 'CellGeneCount', 'ClassName', 'Prob']
                           )
    iss_df.set_index(['Cell_Num'])

    # Ignore the first row. It is the pseudocell to keep the misreads (ie the background)
    iss_df = iss_df[1:]
    summary_logger.info('Data collected!')
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

    out = pd.DataFrame({'Gene': spots.data.gene_name.tolist(),
                        'Gene_id': spots.gene_id.tolist(),
                        'x': spots.data.x.tolist(),
                        'y': spots.data.y.tolist(),
                        'neighbour': max_nbrs,
                        'neighbour_array': nbrs,
                        'neighbour_prob': p})

    return out


def collect_data(cells, spots, genes, single_cell):
    '''
    Collects data for the viewer
    :param cells:
    :param spots:
    :return:
    '''
    iss_df = _iss_summary(cells, genes, single_cell)
    gene_df = _summary(spots)
    return iss_df, gene_df
