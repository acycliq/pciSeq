import pandas as pd
import numpy as np
from distutils import dir_util
from pytest import fixture
import os
from configparser import ConfigParser
from scipy.sparse import load_npz
from pciSeq.src.cell_call import utils
from pciSeq.app import cell_type, fit

@fixture
def datadir(tmpdir, request):
    '''
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    '''
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))

    return tmpdir


def get_scRNAseq(cfg):
    sc_file = cfg.get('dev', 'scRNAseq')

    scRNAseq = pd.read_csv(sc_file, header=None, index_col=0, compression='gzip', dtype=object)
    scRNAseq = scRNAseq.rename(columns=scRNAseq.iloc[0], copy=False).iloc[1:]
    scRNAseq = scRNAseq.astype(np.float).astype(np.uint32)
    return scRNAseq


def get_expected_spots(cfg):
    spots_file = cfg.get('dev', 'spots_file')
    return pd.read_csv(spots_file)


def get_expected_label_image(cfg):
    coo_file = cfg.get('dev', 'coo_file')
    coo_file = utils.load_from_url(coo_file)
    return load_npz(coo_file)


def get_expected_cellData(cfg):
    cellData_file = cfg.get('dev', 'cellData')
    df = pd.read_csv(cellData_file, index_col=0)
    return clean_cellData(df.copy())


def clean_cellData(df):
    tmp_1 = []
    tmp_2 = []
    tmp_3 = []
    for index, row in df.iterrows():
        tmp_1.append(eval(row['Genenames']))
        tmp_2.append(eval(row['CellGeneCount']))
        tmp_3.append(eval(row['ClassName']))
    df['Genenames'] = tmp_1
    df['CellGeneCount'] = tmp_2
    df['ClassName'] = tmp_3
    return df


def get_expected_geneData(cfg):
    geneData_file = cfg.get('dev', 'geneData')
    df = pd.read_csv(geneData_file, index_col=0)
    return clean_geneData(df.copy())


def clean_geneData(df):
    tmp_1 = []
    tmp_2 = []
    for index, row in df.iterrows():
        tmp_1.append(eval(row['neighbour_array']))
        tmp_2.append(eval(row['neighbour_prob']))
    df['neighbour_array'] = tmp_1
    df['neighbour_prob'] = tmp_2
    return df


def test_pciSeq(datadir):
    cfg = ConfigParser()
    cfg.read(datadir.join('config.ini'))
    scRNAseq = get_scRNAseq(cfg)

    iss_spots = get_expected_spots(cfg)
    coo = get_expected_label_image(cfg)

    expected_cellData = get_expected_cellData(cfg)
    expected_geneData = get_expected_geneData(cfg)

    cellData, geneData = fit(iss_spots, coo, scRNAseq)

    a = expected_cellData[['Cell_Num', 'X', 'Y']]
    b = cellData[['Cell_Num', 'X', 'Y']]
    tol = 1e-10
    assert np.max(abs(a.values - b.values)) < tol

    a = cellData[['CellGeneCount', 'Genenames', 'ClassName']]
    b = expected_cellData[['CellGeneCount', 'Genenames', 'ClassName']]
    assert a.equals(b)

    assert expected_geneData.equals(geneData)


