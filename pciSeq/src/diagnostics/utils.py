import os
import redis
import pickle
from ctypes import windll
from pciSeq.src.cell_call.log_config import logger
import pciSeq.src.diagnostics.config as config
import subprocess
import sysconfig


class redis_db():
    def __init__(self, flush=False):
        self.pool = None
        self.redis_client = None
        self._attach_client()
        if flush:
            self._flushdb()

    def _attach_client(self):
        self.pool = config.SETTINGS['POOL']
        self.redis_client = redis.Redis(connection_pool=self.pool)
        assert self.is_connected()

    def to_redis(self, df_in, key, **kwargs):
        df = df_in.copy()
        df = df.assign(**kwargs)
        self.redis_client.set(key, pickle.dumps(df))

    def from_redis(self, key):
        """Retrieve Numpy array from Redis key 'key'"""
        return pickle.loads(self.redis_client.get(key))

    def get_db_tables(self):
        bytestings = self.redis_client.keys()
        return [d.decode('UTF-8') for d in bytestings]

    def _flushdb(self):
        self.redis_client.flushdb()

    def is_connected(self):
        print(os.getcwd())
        try:
            return self.redis_client.ping()
        except (redis.exceptions.ConnectionError, ConnectionRefusedError):
            logger.info("Redis ping failed!. Trying to install redis server")
            os.path.join(os.path.dirname(os.path.realpath(__file__)), )
            msi = os.path.join(sysconfig.get_path('purelib'), 'pciSeq', 'static', 'memurai', 'Memurai-Developer-v3.1.4.msi')
            if not os.path.isfile(msi):
                msi = os.path.join(os.getcwd(), 'static', 'memurai', 'Memurai-Developer-v3.1.4.msi')
            logger.info("Calling %s" % msi)
            exit_code = subprocess.call(msi, shell=True)  # returns the exit code in unix. if 0 then success
            if exit_code > 0:
                raise Exception("Installation failed with exit code: %d" % exit_code)
            return exit_code == 0






