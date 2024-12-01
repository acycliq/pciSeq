import pytest
import pathlib
import hashlib
import numpy as np
from pciSeq.src.core.main import VarBayes
from pciSeq.app import parse_args, stage_data
from pciSeq.src.validation.config_manager import ConfigManager
from pciSeq.src.validation.input_validation import InputValidator
from constants import EXPECTED_AREA_METRICS, EXPECTED_ITER_DELTAS
from utils import setup_logger

test_app_logger = setup_logger(__name__)


def calculate_checksum(str_path):
    """
    Calculate SHA256 hash/checksum of a file
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


class TestPciSeq:
    """Test suite for data processing functionality"""

    def test_parse_args(self, read_demo_data):
        """Test argument parsing with various input combinations"""
        test_app_logger.info('Testing argument parsing')
        spots = read_demo_data[0]
        coo = read_demo_data[1]

        with pytest.raises(ValueError) as excinfo:
            parse_args(spots)
        assert str(excinfo.value) == ('Need to provide the spots and the coo matrix either as keyword arguments or as the first and second positional arguments.')

        _, _, scData, opts = parse_args(spots, coo)
        assert scData is None
        assert opts is None

    def test_validate(self, read_demo_data):
        """Test data validation functionality"""
        test_app_logger.info('Testing data validation')
        spots = read_demo_data[0]
        coo = read_demo_data[1]
        scData = read_demo_data[2]

        with pytest.raises(TypeError) as excinfo:
            cfg_man = ConfigManager.from_opts(None)
            InputValidator.validate(coo, coo, scData, cfg_man)
        assert str(excinfo.value) == "Spots should be passed-in as a dataframe"

    @pytest.mark.parametrize('filename, expected', [
        ('read_demo_data', EXPECTED_AREA_METRICS)
    ])
    def test_stage_data(self, filename, expected, request):
        """Test data staging and processing"""
        test_app_logger.info('Testing data staging and preprocessing')
        read_demo_data = request.getfixturevalue(filename)
        spots = read_demo_data[0]
        coo = read_demo_data[1]

        cells, cell_boundaries, spots = stage_data(spots, coo)
        pytest.fspots = spots
        pytest.fcells = cells

        assert len(cells.label) == len(np.unique(cells.label))
        assert cells.area.sum() == expected['area_sum']
        assert cells.label.max() == expected['max_label']
        assert np.all(
            spots[['x_global', 'y_global', 'label', 'x_cell', 'y_cell']]
            .sum().round(5).values == expected['column_sums'])

    @pytest.mark.parametrize('filename, expected', [
        ('read_demo_data', EXPECTED_ITER_DELTAS)
    ])
    def test_varBayes(self, filename, expected, request):
        """Test VarBayes algorithm convergence"""
        test_app_logger.info('Testing algorithm convergence')
        read_demo_data = request.getfixturevalue(filename)
        spots, coo, scData = read_demo_data

        # Create config
        opts = {
            'launch_viewer': True,
            'launch_diagnostics': True,
            'max_iter': 31,
        }
        cfg_man = ConfigManager.from_opts(opts)

        # validate inputs
        spots, coo, scdata, cfg = InputValidator.validate(spots, coo, scData, cfg_man)

        _cells, cellBoundaries, _spots = stage_data(spots, coo)
        varBayes = VarBayes(_cells, _spots, scData, cfg)
        cellData, geneData = varBayes.run()

        arr_1 = np.array(varBayes.iter_delta, dtype=np.float32).round(11)
        arr_2 = np.array(expected[: cfg['max_iter']], dtype=np.float32).round(11)
        assert np.allclose(arr_1, arr_2, rtol=0, atol=1e-05)
