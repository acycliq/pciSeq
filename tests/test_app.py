import os
import numpy as np
import pathlib
import hashlib
import pytest
from pciSeq.src.core.utils import get_out_dir
from pciSeq.src.core.main import VarBayes
import pciSeq.config as config
from pciSeq.app import parse_args, validate, stage_data
import logging


# set up the logger
# logger_setup()


def calculate_checksum(str_path):
    """
    Calculate SHA256 hash/checksum of a file

    Args:
        path (str): Path to file

    Returns:
        str: checksum as hexadecimal string

    From :
    https://github.com/metoppv/improver/blob/master/improver_tests/acceptance/acceptance.py
    """
    hasher = hashlib.sha256()
    path = pathlib.Path(str_path)
    with open(path, mode="rb") as kgo_file:
        while True:
            # read 1 megabyte binary chunks from file and feed them to hasher
            kgo_chunk = kgo_file.read(2 ** 20)
            if not kgo_chunk:
                break
            hasher.update(kgo_chunk)
    checksum = hasher.hexdigest()
    return checksum


def tmpdir():
    out_dir = os.path.join(get_out_dir(), 'tests')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_dir


def test_parse_args(read_demo_data):
    logging.getLogger().info('test_2')
    spots = read_demo_data[0]
    coo = read_demo_data[1]

    with pytest.raises(AssertionError) as excinfo:
        parse_args(spots)
    assert str(excinfo.value) == ('Need to provide the spots and the coo matrix as the first '
                                  'and second args to the fit() method.')

    _, _, scData, opts = parse_args(spots, coo)
    assert scData is None
    assert opts is None


def test_validate(read_demo_data):
    logging.getLogger().info('test_2')
    spots = read_demo_data[0]
    coo = read_demo_data[1]

    with pytest.raises(AssertionError) as excinfo:
        validate(coo, coo, config.DEFAULT)
    assert str(excinfo.value) == ("Spots should be passed-in to the fit() method as "
                                  "a dataframe with columns ['Gene', 'x', 'y']")


def test_stage_data(read_demo_data):
    spots = read_demo_data[0]
    coo = read_demo_data[1]
    cells, cell_boundaries, spots = stage_data(spots, coo)
    pytest.fspots = spots
    pytest.fcells = cells
    assert len(cells.label) == len(np.unique(cells.label))
    assert cells.area.sum() == 3721912.0
    assert cells.label.max() == 3481
    assert np.all(spots[['x_global', 'y_global', 'label', 'x_cell', 'y_cell']].sum().round(5).values == [231498829.00000, 150952363.00000, 75868924.00000, 151411235.74292, 82555694.21166])


def test_varBayes(read_demo_data):
    scData = read_demo_data[2]
    varBayes = VarBayes(pytest.fcells, pytest.fspots, scData, config.DEFAULT)
    varBayes.run()

    expected_delta = np.array([1.0, 0.8998741535450602, 0.6707994416612973, 0.4703557918803584, 0.6274678470436773, 0.5809977507020871,
     0.6050537780421158, 0.281664156778341, 0.26291180844419526, 0.1911332823201044, 0.5582350505154667,
     0.1904928001599277, 0.2482987581795233, 0.2825277025986648, 0.06827397193459617, 0.16097980010421764,
     0.5828276333586281, 0.16734418317139443, 0.14889440752498986, 0.1229702319182493, 0.08985526872588012,
     0.044863509416499414, 0.0678009601806856, 0.10751329267974208, 0.19152470052524656, 0.16275473631515358,
     0.06207517767135862, 0.04486396122645342, 0.0315151843141393, 0.01963733267979928]).round(11)

    assert np.all(np.array(varBayes.iter_delta).round(11) == expected_delta)



# # def read_dataframe():
# #     """
# #     This is for the test at the bottom, needs a lot more work
# #     """
# #     cellData_path = os.path.join(get_out_dir(), 'pciSeq', 'data', 'cellData.tsv')
# #     cellData = pd.read_csv(cellData_path, sep='\t')
# #     cellData[['Cell_Num', 'X', 'Y']].round(6)
# #     cellData['CellGeneCount'] = cellData['CellGeneCount'].apply(lambda x: np.round(np.array(eval(x)), 3))
# #     cellData['Prob'] = cellData['Prob'].apply(lambda x: np.round(np.array(eval(x)), 3))
# #     return cellData
#
#
# bbox = [
#     (4238, 364),  # bottomleft, [x0, y0]
#     (5160, 933)  # topright, [x1, y1]
# ]
#
# spots, image_label = make_test_data(tmpdir(), bbox=None)
# coo = coo_matrix(image_label)
#
# scRNAseq = pd.read_csv('test_scRNAseq.csv').set_index('gene_name')
#
# # main task
# # _opts = {'max_iter': 10}
# opts = {'save_data': True,
#          'launch_viewer': False,
#          'launch_diagnostics': False,
#          'output_path': [tmpdir()]
#         }
#
# test_list = [
#     spots,
#     coo,
#     scRNAseq,
#     opts
# ]
#
# expected_list = ['364e544d2c837902645083b8c2b9298fef2987a24c85160ec61c50d2e6f7ffce']
#
#
# @pytest.mark.parametrize('test_data, expected', [
#     (test_list, expected_list)
# ])
# def test_app(test_data, expected):
#     expected_csum = expected[0]
#     cellData, geneData = fit(spots=spots, coo=coo, scRNAseq=scRNAseq, opts=opts)
#
#     csum_pickle = calculate_checksum(os.path.join(pathlib.Path(tmpdir()), 'pciSeq', 'data', 'debug', 'pciSeq.pickle'))
#     print(csum_pickle)
#     assert csum_pickle == expected_csum
#
#     # Clean the temporary dir
#     shutil.rmtree(tmpdir())




# THE TEST BELOW NEEDS MORE WORK
# test_cellData_df = read_dataframe()
# expected_cellData_df = read_dataframe()


# @pytest.mark.parametrize('test_cellData, expected_cellData', [
#     (test_cellData_df, expected_cellData_df)
# ])
# def test_celldata_df(test_cellData, expected_cellData):
#     ## This needs more work. it is not finished!
#     assert assert_frame_equal(test_cellData, expected_cellData)