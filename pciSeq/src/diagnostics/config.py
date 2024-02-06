import redis

SETTINGS = {
    # 'DB_URL': os.path.join(get_out_dir(pciSeq_cfg._BASE['output_path'], ''), 'pciSeq.db'),
    'POOL': redis.ConnectionPool(host='localhost', port=6379, db=0)
}

