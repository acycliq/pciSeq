import numpy as np
import pandas as pd
import os
import xarray as xr
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
CONFIG_FILE = dir_path + '/config.yml'


logger = logging.getLogger()


def _load_geneset(config):
    '''
    :param config:
    :return:
    '''

    sc_path = os.path.join(dir_path, config['scRNAseq'])
    logger.info('Loading geneset from %s' % sc_path)
    df = pd.read_csv(sc_path, header=None, index_col=0, compression='gzip')

    # from https://github.com/pandas-dev/pandas/issues/19383
    df = df.rename(columns=df.iloc[0], copy=False).iloc[1:]

    # make a dataframe from the xarray
    logger.info('Renaming subclasses PC.CA2 and PC.CA3 to be PC.Other1 and PC.Other2')
    class_name = ['PC.Other1' if x == 'PC.CA2' else x for x in df.columns.values]
    class_name = ['PC.Other2' if x == 'PC.CA3' else x for x in class_name]

    df.columns = class_name

    df = df.rename_axis("class_name", axis="columns") \
        .rename_axis('gene_name') \
        .astype('float64')

    return xr.DataArray(df)


def _scloader(config):
    sc_path = os.path.join(dir_path, config['scRNAseq'])
    df = pd.read_csv(sc_path, header=0, index_col=0, compression='gzip')\
        .rename_axis("Class", axis="columns")\
        .rename_axis('GeneName')\
        .astype('float64')
    return df



def _normalise(df):
    '''
    removes columns with zero mean (ie all zeros) and then rescales so that
    total counts remain the same.
    :param df:
    :return:
    '''

    # find the mean for each column
    col_mean = df.mean(axis=0)

    isZero = col_mean == 0.0

    # remove column if its mean is zero
    df2 = df.loc[:, ~isZero]

    total = df.sum().sum()
    total2 = df2.sum().sum()

    # rescale df2 so that the total counts are the same
    df2 = df2 * total/total2

    # sanity check
    assert df.sum().sum() == df2.sum().sum()

    return df2


def _remove_zero_cells(da):
    '''
    Effectivelly what this does is to simply remove zero columns
    :param da:
    :return:
    '''
    da = da.loc[:, (da != 0).any('gene_name')]

    return da


def geneSet(spots, config):
    genes = spots.gene_panel.index.values
    da = _load_geneset(config)

    # I should be able to just do: da.loc[genes] but it gives an error
    # and I dont know why!!
    # Hence i am doing this funny workaround
    da = xr.DataArray(da.to_pandas()
                        .loc[genes]
                        .rename_axis("class_name", axis="columns")
                        .rename_axis('gene_name')
                        )

    # da = _normalise(da)
    da = _remove_zero_cells(da)

    # calc the mean expression within each cell type
    mean_expression_da = config['Inefficiency'] * da.groupby('class_name').mean(dim='class_name')

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
    log_mean_expression = np.log(mean_expression + config['SpotReg'])

    ds = xr.Dataset({'mean_expression': mean_expression,
                      'log_mean_expression': log_mean_expression
                     })

    return ds
