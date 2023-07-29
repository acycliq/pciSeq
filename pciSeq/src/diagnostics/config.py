from pciSeq.src.core.utils import get_out_dir
import pciSeq.config as pciSeq_cfg
import redis
import os


SETTINGS = {
    # 'DB_URL': os.path.join(get_out_dir(pciSeq_cfg._BASE['output_path'], ''), 'pciSeq.db'),
    'POOL': redis.ConnectionPool(host='localhost', port=6379, db=0)
}

