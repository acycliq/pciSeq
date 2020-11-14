import numpy as np
from src.preprocess.utils import load_mat, split_CellMap
import urllib.request
import base64
import json
from scipy.sparse import coo_matrix
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


class Tile:
    def __init__(self, cfg):
        self.cfg = cfg
        self._cellmap_chunks = split_CellMap(cfg['cellmap_full'], cfg['tile_size'], cfg['tile_size'])
        self.tiles_across = -(-cfg['img_width'] // cfg['tile_size'])
        self.tiles_down = -(-cfg['img_height'] // cfg['tile_size'])
        self.tile_shape = [cfg['tile_size'], cfg['tile_size']]
        self.tiles = self.populate_tiles(self.tiles_across * self.tiles_down)
        del self._cellmap_chunks

    def populate_tiles(self, tiles_counts):
        """
        Populates the tiles. Each fov is a dict with key/value pairs:
            fov_id: (int) The id of the fov
            fov_range: a list of coordinates [[x_min, x_max], [y_min, y_max]] of the x-side and y-side of the fov
            fov_offset_x: (int) the x-coordinate of the top left corner of the fov
            fov_offset_y: (int) the y-coordinate of the top left corner of the fov

        Parameters
        ----------
        tiles_counts: (int) The counts of tiles in the image

        Returns
        -------
        a list of tiles
        """

        logger.info('loading specs for each tile..')
        n = max(1, cpu_count() - 1)
        pool = ThreadPool(n)
        results = pool.map(self._helper, range(tiles_counts))
        pool.close()
        pool.join()
        logger.info('Tile specs loaded!')

        tile_attr = []
        for i, d in enumerate(results):
            assert i == d[0]['tile_id'], 'data are not aligned'
            tile_attr.append(d[0])

        shape_x = np.unique([d[1] for d in results])
        shape_y = np.unique([d[2] for d in results])

        assert len(shape_x) == 1, 'All tiles should have the same x-length'
        assert len(shape_y) == 1, 'All tiles should have the same y-length'
        assert [shape_x[0], shape_y[0]] == self.tile_shape, 'Shape is wrong'
        return tile_attr

    def _helper(self, f):
        coords = self.get_tile_coords(f)
        temp = {'tile_id': f,
                'tile_range': {'x': coords[0], 'y': coords[1]},
                'tile_offset_x': coords[0][0],
                'tile_offset_y': coords[1][0],
                'label_image': self.get_label_img(f)
                # '_label_image_outlines': None,
                }
        return temp, coords[0][1] - coords[0][0], coords[1][1] - coords[1][0]

    def get_tile_origin(self, fov_id):
        """
        return the coordinates of the top left corner of a tile given its id

        Parameters
        ----------
        fov_id: ths id of the fov

        Returns
        -------
        the coordinates of the top left corner
        """
        tiles_across = self.tiles_across
        x_size = self.tile_shape[0]
        y_size = self.tile_shape[1]

        # find how far right you are:
        x = fov_id % tiles_across

        # find how far down you are:
        y = fov_id // tiles_across
        x0 = x_size * x
        y0 = y_size * y
        return x0, y0

    def get_tile_coords(self, fov_id):
        """
        return the range of the x and y sides of the fov given the fov id.
        The ranges are in global coordinates

        :param fov_id:
        :return:
        """
        x, y = self.get_tile_origin(fov_id)
        x_size = self.tile_shape[0]
        y_size = self.tile_shape[1]

        x_range = [x, x + x_size]
        y_range = [y, y + y_size]

        return x_range, y_range

    # def npz_path(self, i, ini):
    #     return os.path.join(ini['cellpose_out'], 'fov_%d.tif.npz' % i)

    def get_label_img(self, i):
        return coo_matrix(self._cellmap_chunks[i])
        # # npz_file = self.npz_path(i, ini)
        # # b = np.load(npz_file)
        # cellmap = load_mat(filepath)
        # coo = coo_matrix((b['data'], (b['row'], b['col'])), shape=b['shape'])
        # return coo


if __name__ == "__main__":
    fov = Tile(config.ini)

    print('Done!')