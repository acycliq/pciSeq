from pciSeq.src.cell_call.utils import get_out_dir
import pciSeq.config as pciSeq_cfg
import os


SETTINGS = {
    # 'DB_URL': "file:memdb1?mode=memory&cache=shared"
    # 'DB_URL': r'D:\Home\Dimitris\OneDrive - University College London\dev\Python\pciSeq\pciSeq\my_pciSeq.db',
    'DB_URL': os.path.join(get_out_dir(pciSeq_cfg.DEFAULT['output_path'], ''), 'pciSeq.db')
}

