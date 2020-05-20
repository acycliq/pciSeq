import numpy as np
import urllib.request
import base64
import json
import os
from multiprocessing.dummy import Pool as ThreadPool
import config
import logging
import time

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

logger = logging.getLogger()
# logger.disabled = True

class Fov:
    def __init__(self, fovs_across, fovs_down, cfg):
        self.fovs_across = fovs_across
        self.fovs_down = fovs_down
        self.my_config = cfg
        self.fovs = []
        self.scaling_factor = 1
        self.fov_shape = None
        self.pool = None
        start = time.time()
        self.fov_dict_par(fovs_across * fovs_down)
        print('Finished parallel in %s' % (time.time() - start))

    def setThreadPool(self, n):
        self.pool = ThreadPool(n)

    def fov_dict_par(self, fovs_num):
        logger.info('loading specs for each fov..')
        pool = ThreadPool(14)
        results = pool.map(self._helper, range(fovs_num))
        pool.close()
        pool.join()
        logger.info('Fov specs loaded!')

        self.fovs = [d[0] if d[0]['fov_id']==i else 1/0 for i, d in enumerate(results)]
        _fov_shape_x = [d[1] for d in results]
        _fov_shape_y = [d[2] for d in results]

        assert len(set(_fov_shape_x)) == 1
        assert len(set(_fov_shape_y)) == 1
        self.fov_shape = np.array([_fov_shape_x[0], _fov_shape_y[0]]).astype(np.int32)

    def _helper(self, f):
        coords = self.get_fov_coords(f)
        # coords = self.global_coords(coords)
        temp = {'fov_id': f,
                'fov_range': {'x': coords[0], 'y': coords[1]},
                'fov_offset_x': coords[0][0],
                'fov_offset_y': coords[1][0],
                'fov_scaling_factor': 1,
                'fov_dir': self.get_dir(f),
                '_label_image_outlines': None,
                }
        return temp, coords[0][1] - coords[0][0], coords[1][1] - coords[1][0]

    def global_coords(self, coords):
        x_range = coords[0]
        y_range = coords[1]
        k = coords[2]
        x_global = [d * k for d in x_range]
        y_global = [d * k for d in y_range]
        return x_global, y_global, k

    def get_fov_origin(self, fov_id):
        '''
        return the coordinates of the top left corner of a tile given its id

        :param fov_id:
        :return:
        '''
        fovs_across = self.fovs_across
        fov_size = self.my_config['fov_size']

        # find how far right you are:
        x = fov_id % fovs_across

        # find how far down you are:
        y = fov_id // fovs_across

        x0 = fov_size * x
        y0 = fov_size * y

        return x0, y0


    def get_fov_coords(self, fov_id):
        '''
        return the range of the x and y sides of the fov given the fov id.
        The ranges are in global coordinates

        :param fov_id:
        :return:
        '''
        x, y = self.get_fov_origin(fov_id)
        fov_size = self.my_config['fov_size']

        x_range = [x, x + fov_size]
        y_range = [y, y + fov_size]

        return x_range, y_range


    def get_dir(self, fov_id):
        return os.path.join(self.my_config['FOV_ROOT'], 'fov_' + str(fov_id))


if __name__ == "__main__":
    fovs_across = 11
    fovs_down = 14

    fov = Fov(fovs_across, fovs_down, config.PREPROCESSOR)

    print('Done!')
