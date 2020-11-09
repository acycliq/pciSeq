import numpy as np
import pandas as pd
import skimage.measure as skmeas
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)
logger = logging.getLogger()


class Post_merge:
    def __init__(self, tiles, merge_register, label_generator):
        self.tiles = tiles
        self.merge_register = merge_register
        self.label_generator = label_generator

    def run(self):
        """ Assign new labels to the cells. I will also create a new key '''image_objects''' in the fov dict
        with value a dataframe that keeps the label, fov_id, area, and centroid (local) coords of the cells
        """
        self.process_tiles()
        # Need also to flip the sign of the labels in the merge_register dict
        self.flip_sign()

    def process_tiles(self):
        n = max(1, cpu_count() - 1)
        pool = ThreadPool(n)
        results = pool.map(self.process_tile, self.tiles)
        pool.close()
        pool.join()
        return results

    def process_tile(self, tile):
        # set the global labels
        self.global_labels(tile)

        # Finally, get area, centroid coords and label of the cells you find on the fov
        self.image_objects(tile)
        return True

    def global_labels(self, tile):
        # logger.info('tile_id is: %d' % tile['tile_id'])
        data = tile['label_image'].data

        # Get the positive labels. Negative means the cell is clipped, hence it has been
        # given a new label (negative valued)
        labels = np.unique(data[data > 0])
        label_map = {d: self.label_generator() for d in labels}
        tile['label_image'].data = np.array([label_map[d] if d > 0 else d for d in data])

        # Sanity check. merge register keeps the cells that are clipped, the cells that span
        # over more than a single fov
        # Get the labels from the merge register for this particular fov
        clipped_cell_labels = {k: v for k, v in self.merge_register.entries.items() if tile['tile_id'] in v}.keys()
        assert np.all([d in tile['label_image'].data for d in clipped_cell_labels]), \
            "A label that the cell register says should exist in the fov %d doesnt seem to be verified by the " \
            "label_image" % tile['tile_id']

        # Flip now the sign of the labels
        tile['label_image'].data = -1 * tile['label_image'].data
        if np.any(tile['label_image'].data):
            logger.info('fov:%d, Global labels are set.' % tile['tile_id'])
        else:
            logger.info('fov:%d empty, nothing to do here' % tile['tile_id'])


    def image_objects(self, tile):
        """
        Reads a label_image and returns a dataframe with the objects found on it.
        The array has index the label and columns: fov_id, area and the cell centroids x, y
        """
        tile_id = tile['tile_id']
        label_image = tile['label_image'].toarray().astype(np.int32)

        props = skmeas.regionprops(label_image)
        props_arr = [[int(p.label), int(tile_id), int(p.area), p.centroid[1], p.centroid[0]] for i, p in
                     enumerate(props)]
        props_df = pd.DataFrame(props_arr, columns=['label', 'tile_id', 'area', 'x_local', 'y_local'])

        # stick the df to the dict
        tile['image_objects'] = props_df  # <--- Not sure If I should do this or just return the df and add the dictionary filed at the calliing function
        return tile

    def flip_sign(self):
        # need to flip the sign of the keys in the merge_register dict
        # Maybe I should add some sanity checks here.
        # For example I can do the flipping inside the loop of the calling function and while
        # you loop get the keys (ie the labels) for the particular fov and then flip the label
        # sign (if it negative because it may have been flipped already if the label exists in
        # more that one fov

        if bool(self.merge_register):
            my_keys = list(self.merge_register.entries.keys())
            vals = list(self.merge_register.entries.values())

            assert max(my_keys) < 0, 'All keys should be non-zero negative'
            new_keys = [-1 * d for d in my_keys]

            self.merge_register.entries = dict(zip(new_keys, vals))
        else:
            logger.info('Merge register empty')

