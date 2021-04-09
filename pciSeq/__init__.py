from pciSeq.app import fit
from pciSeq.app import cell_type
from pciSeq.src.preprocess.spot_labels import stage_data
import pciSeq.src.cell_call.utils as utils
import json
from os.path import dirname

with open(dirname(__file__) + '/version.json') as fp:
    _info = json.load(fp)

__version__ = _info['version']