import numpy as np
import pandas as pd
import os
import xarray as xr
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger()


def load_scRNAseq(config):
    """
    :param config:
    :return:
    """

    logger.info('Loading single cell data from %s' % config['scRNAseq'])
    # Exression table from single cell data can have counts from the same class multiple times (more than once) .
    # Pandas do not like headers with non-unique names, hence do not use a header
    # when you read the csv.
    # Read the csv with no headers  The class name will be the very first row of the dataframe
    df = pd.read_csv(config['scRNAseq'], header=None, index_col=0, compression='gzip', dtype=object)

    # from https://github.com/pandas-dev/pandas/issues/19383
    # Move now the class name from the top row to the headers of the dataframe
    df = df.rename(columns=df.iloc[0], copy=False).iloc[1:]

    # The dataframe has now gene name as row labels (index) and class names as headers
    # The values however are still strings. Convert them to uint32
    df = df.astype(np.float).astype(np.uint32)

    # Some housekeeping. PC.CA2 and PC.CA3 should be renamed if they are in the scRNAseq data
    logger.info('Renaming subclasses PC.CA2 and PC.CA3 to be PC.Other1 and PC.Other2')
    class_name = ['PC.Other1' if x == 'PC.CA2' else x for x in df.columns.values]
    class_name = ['PC.Other2' if x == 'PC.CA3' else x for x in class_name]

    df.columns = class_name
    df = df.rename_axis("class_name", axis="columns") \
        .rename_axis('gene_name')
    return df


def _remove_zero_cells(df):
    """
    Removes zero columns (ie if a column is populated by zeros only, then it is removed)
    :param da:
    :return:
    """
    out = df.loc[:, (df != 0).any(axis=0)]
    return out


def sc_expression_data(genes, config):
    gene_names = genes.gene_names
    df = load_scRNAseq(config)
    df = df.loc[gene_names].rename_axis('gene_name')

    df = _remove_zero_cells(df.copy())
    da = xr.DataArray(df)

    # calc the mean expression within each cell type
    mean_expression_da = config.getfloat('Inefficiency') * da.groupby('class_name').mean(dim='class_name')

    # sort the dataframe (the index only)
    mean_expression = mean_expression_da.sortby('gene_name')

    # append the zero cell
    zero_df = pd.DataFrame({'Zero': np.zeros(mean_expression.shape[0])
                            }, index=da.gene_name.values
                           )\
        .rename_axis("class_name", axis="columns")\
        .rename_axis('gene_name')

    zero_da = xr.DataArray(zero_df)
    mean_expression = xr.concat([mean_expression_da, zero_da], 'class_name')
    log_mean_expression = np.log(mean_expression + config.getfloat('SpotReg'))

    ds = xr.Dataset({'mean_expression': mean_expression,
                      'log_mean_expression': log_mean_expression
                     })

    return ds
