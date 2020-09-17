import numpy as np
import urllib.request
import base64
import json
import os
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
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
    '''
    The Fov object contains the fovs

    Parameters
    ----------
    cfg: dict
        dictionary with key/value pairs:
            fovs_across: (int) How many fovs the full image has from left to right
            fovs_down: (int) How many fovs the full image has from top to bottom
            fov_shape: (list) list-like of length 2 describing the size in px of the fov.
                        It is assumed that all fovs have the same shape. If fov_shape = [1000, 2000]
                        then it is implied that the x side of the fov is 1000px long and the y-side 2000px
    '''
    def __init__(self, cfg):
        self.fovs_across = cfg['fovs_across']
        self.fovs_down = cfg['fovs_down']
        self.fov_shape = cfg['fov_shape']
        self.fovs = self.populate_fovs(self.fovs_across * self.fovs_down)
        self.scaling_factor = 1  # This will be removed
        start = time.time()
        print('Finished parallel in %s' % (time.time() - start))

    def populate_fovs(self, fovs_counts):
        """
        Populates the fovs. Each fov is a dict with key/value pairs:
            fov_id: (int) The id of the fov
            fov_range: a list of coordinates [[x_min, x_max], [y_min, y_max]] of the x-side and y-side of the fov
            fov_offset_x: (int) the x-coordinate of the top left corner of the fov
            fov_offset_y: (int) the y-coordinate of the top left corner of the fov

        Parameters
        ----------
        fovs_counts: (int) The counts of fovs in the image

        Returns
        -------
        a list of fovs
        """

        logger.info('loading specs for each fov..')
        n = max(1, cpu_count() - 1)
        pool = ThreadPool(n)
        results = pool.map(self._helper, range(fovs_counts))
        pool.close()
        pool.join()
        logger.info('Fov specs loaded!')

        fov_attr = []
        for i, d in enumerate(results):
            assert i == d[0]['fov_id'], 'data are not aligned'
            fov_attr.append(d[0])

        shape_x = np.unique([d[1] for d in results])
        shape_y = np.unique([d[2] for d in results])

        assert len(shape_x) == 1, 'All fovs should have the same x-length'
        assert len(shape_y) == 1, 'All fovs should have the same y-length'
        assert [shape_x[0], shape_y[0]] == self.fov_shape, 'Shape is wrong'
        return fov_attr

    def _helper(self, f):
        coords = self.get_fov_coords(f)
        temp = {'fov_id': f,
                'fov_range': {'x': coords[0], 'y': coords[1]},
                'fov_offset_x': coords[0][0],
                'fov_offset_y': coords[1][0],
                'fov_scaling_factor': 1,
                # '_label_image_outlines': None,
                }
        return temp, coords[0][1] - coords[0][0], coords[1][1] - coords[1][0]

    def get_fov_origin(self, fov_id):
        """
        return the coordinates of the top left corner of a tile given its id

        Parameters
        ----------
        fov_id: ths id of the fov

        Returns
        -------
        the coordinates of the top left corner
        """
        fovs_across = self.fovs_across
        x_size = self.fov_shape[0]
        y_size = self.fov_shape[1]

        # find how far right you are:
        x = fov_id % fovs_across

        # find how far down you are:
        y = fov_id // fovs_across
        x0 = x_size * x
        y0 = y_size * y
        return x0, y0

    def get_fov_coords(self, fov_id):
        """
        return the range of the x and y sides of the fov given the fov id.
        The ranges are in global coordinates

        :param fov_id:
        :return:
        """
        x, y = self.get_fov_origin(fov_id)
        x_size = self.fov_shape[0]
        y_size = self.fov_shape[1]

        x_range = [x, x + x_size]
        y_range = [y, y + y_size]

        return x_range, y_range


if __name__ == "__main__":
    fovs_across = 11
    fovs_down = 14

    fov = Fov(fovs_across, fovs_down, config.PREPROCESSOR)

    print('Done!')
