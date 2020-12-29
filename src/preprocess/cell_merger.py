import os
import shutil
import logging
import itertools
import numpy as np
from collections import defaultdict
from src.preprocess.post import Post_merge
from src.preprocess.utils import _to_csr_matrix, _get_connected_labels
from scipy.sparse.csgraph import connected_components

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)
logger = logging.getLogger()


class Merge_register(object):
    def __init__(self, parent):
        self.entries = defaultdict(list)
        self.parent = parent

    def update_register(self, tile_id, label, old_label):
        self.entries[label].append(tile_id)
        self.entries[label] = sorted(list(set(self.entries[label])))

        logger.info('tile_%d: label %d ---> label %d' % (tile_id, old_label, label))
        if (old_label != label) and (old_label < 0):
            self.replace_label(old_label, label)

    def replace_label(self, old_label, label):
        # replace the register
        _dict = self.entries
        for tile_id in _dict[old_label]:
            if tile_id not in _dict[label]:
                _dict[label].append(tile_id)
            mask = self.parent.tiles[tile_id]['label_image'].data == old_label
            self.parent.tiles[tile_id]['label_image'].data[mask] = label
            logger.info('tile: %d: replaced labels in "label_image" that were equal to %d with %d' % (tile_id, old_label, label))

        logger.info('Dropped key, value pair: (%d, %s) from the merge_register' % (old_label, _dict[old_label]))
        _dict.pop(old_label)

# ----------------------------------------------------------------------------------------------------------------------

class Stage(object):
    def __init__(self, tile_obj, spots_all):
        self.counter = itertools.count()
        self.merge_register = Merge_register(self)
        self.cell_props = None
        self.spots = None
        self.tiles = tile_obj.tiles
        self.tiles_across = tile_obj.tiles_across
        self.tiles_down = tile_obj.tiles_down
        self.tile_shape = tile_obj.tile_shape
        self.scaling_factor = 1  # ti
        self.compose_dict(spots_all)

    def compose_dict(self, spots_all):
        """
        Mutates in-place the tile object by adding two more key/value pairs

        Parameters
        ----------
        spots_all: dataframe
            Contains all the spots (Gene names and x,y coords) for the full image
        """
        for i, d in enumerate(self.tiles):
            # d['label_image'] = self.tile_label_image(i)  # label_image for the i-th tile
            d['spots'] = self.tile_spots(spots_all, i)   # spots for the i-th tile

    def tile_spots(self, data, i):
        """ spots for the i-th tile """
        x_range = self.tiles[i]['tile_range']['x']
        y_range = self.tiles[i]['tile_range']['y']
        mask = (data.x.values >= x_range[0]) & \
               (data.x.values < x_range[1]) & \
               (data.y.values >= y_range[0]) & \
               (data.y.values < y_range[1])
        df = data[mask].dropna()

        df = df[['Gene', 'x', 'y']]
        df = df[~df.duplicated()]
        gene_name, idx = np.unique(df.Gene.values, return_inverse=True)
        df['gene_id'] = idx  # this is effectively the gene id

        df['x'] = df.x * self.scaling_factor
        df['y'] = df.y * self.scaling_factor
        df = df.sort_values(['x', 'y'], ascending=[True, True]) \
            .reset_index(drop=True)  # <-- DO NOT FORGET TO RESET THE INDEX
        return df

    def post_merge(self, argin):
        pm = Post_merge(argin[0], argin[1], argin[2])
        pm.run()

    def label_generator(self):
        return -1 * (next(self.counter) + 1)

    def merge_cells(self):
        """  Merge cells clipped by two or more tiles. """
        for tile in self.tiles:
            logger.info('\n')
            logger.info('Doing tile %i' % tile['tile_id'])
            self.merge(tile)
        logger.info('Relabelling finished')
        logger.info('\n')

    def merge(self, tile):
        """
        Does most of the heavy lifting for cell merging. Mutates in-place the label_image arrays of three tiles.
        If tile has tile_id = i then the mutated label_images are for the tiles with:
            tile_id = i
            tile_id = i + 1 (the neighbouring tile at the right)
            tile_id = i - #tiles_across (the neighbouring tile at the top)
        Parameters
        ----------
        tile: an instance of the class Fov

        Notes
        -----
        Is is assumed that each tile is big enough (relative to the cells) so that there is no cell bigger in size that a tile.
        For example, the most complicated case will be a cell clipped by four tiles forming a 2x2 setup with the cell centroid close
        at the intersection of the four tiles
        """
        tile_id = tile['tile_id']
        adj_img = self.adjacent_tile(tile_id)

        logger.info('tile_%d neighbours: (above, left): (%s, %s)' % (tile_id, adj_img['up'], adj_img['left']))

        # Bottom border of the label array above
        if (adj_img['up'] is not None) and np.any(self.tiles[adj_img['up']]['label_image'].data):
            tile_up = self.tiles[adj_img['up']]
            coo_aa, coo_bb = self.dissolve_borders(tile_up, tile, transpose=True)
            tile_up['label_image'] = coo_aa
            tile['label_image'] = coo_bb

        if adj_img['left'] is not None:
            tile_left = self.tiles[adj_img['left']]
            coo_a, coo_b = self.dissolve_borders(tile_left, tile)
            tile_left['label_image'] = coo_a
            tile['label_image'] = coo_b

    def adjacent_tile(self, tile_id):
        if tile_id % self.tiles_across != 0:
            left = tile_id - 1
        else:
            left = None

        if tile_id >= self.tiles_across:
            up = tile_id - self.tiles_across
        else:
            up = None
        return {'left': left, 'up': up}

    def dissolve_borders(self, adjc_tile, tile, transpose=False):
        """
        Compares the label_image arrays from two neighbouring (one next another) tiles. If the last column of the
        label_image at the left and the first column of the one at the right have non-zero values at the same location
        then the labels at these locations are assigned a new and common label
        Parameters
        ----------
        adjc_tile: an instance of the class Fov
            The neighbouring tile. Could be the neighbour from the right, or from above
        tile: an instance of the class Fov
            the current tile
        transpose: bool. Optional
            if adjc_tile is the neighbour from the top, then set this to True. Default is False


        Returns
        -------
        temp_a, temp_b: tuple
            A tuple of two label_image arrays that correspond to the adjacent and the current tile respectively
        """
        if transpose:
            adjc_img = adjc_tile['label_image'].transpose()
            img = tile['label_image'].transpose()
        else:
            adjc_img = adjc_tile['label_image']
            img = tile['label_image']

        arr = adjc_img.toarray()
        adjc_border = arr[:, -1]
        border = img.toarray()[:, 0]

        logger.info('length of adjc_border: %d' % adjc_border.shape[0])
        logger.info('length of adjc_border: %d' % border.shape[0])
        matched_labels = self.connect_labels(adjc_border, border)
        temp_a = self.tiles[adjc_tile['tile_id']]['label_image'].copy()
        temp_b = self.tiles[tile['tile_id']]['label_image'].copy()

        for d in matched_labels:
            new_label = self._new_label(d)
            for x in d['a']:
                temp_a.data[temp_a.data == x] = new_label
                self.merge_register.update_register(adjc_tile['tile_id'], new_label, x)
                # logger.info('tile_%d: label %d ---> label %d' % (adjc_tile['tile_id'], x, new_label))

            for x in d['b']:
                temp_b.data[temp_b.data == x] = new_label
                self.merge_register.update_register(tile['tile_id'], new_label, x)
                # logger.info('tile_%d: label %d ---> label %d' % (tile['tile_id'], x, new_label))
        return temp_a, temp_b

    def _new_label(self, d):
        # get a list from the dict values
        _list = [x[0] for x in list(d.values())]

        # Find the biggest non-positive value
        m = sorted([el for el in _list if el < 0])
        # m = list(set(m))
        if len(m) > 0:
            # label has already been reassigned. Give that to all merging cells and do not generate a new label.
            out = m[-1]
            # I think m should contain the same elements anyway. If len(set(m)) > 1 then something went wrong??
            logger.info('m is: %s' % m)
        else:
            out = self.label_generator()

        assert out < 0, 'Generated labels should be negative'
        return out

    def connect_labels(self, par_a, par_b):
        '''
        compares two list-like input objects of the same size and returns the elements in ''par_a''
        and ''par_b'' which are non-zero and have the same index position in both inputs
        Example connect_labels([0,0,0,2,2,2,4,7], [0,2,2,2,2,2,2,9]) returns
            [
                {'a': [2, 4], 'b': [2]},
                {'a': [7],    'b': [9]}
             ]
        which means that from the first arg the values 2 and 4 meet (have the same position in the array)
        with the value of 2 from the second arg.
        Also value 7 from the first arg has the same position with value 9 in the second arg. They are the
        last elements in both lists
        :param a: list
        :param b: list
        :return:
        '''

        assert len(par_a) == len(par_b), "inputs to the function should have the same length"

        a, b, lookup_label_a, lookup_label_b = self._shift_labels(par_a, par_b)
        assert len(a) == len(b), "a and b do not have the same length"
        assert len(a) == len(par_a)
        assert len(b) == len(par_b)
        # Make sure the two list do not have common elements
        a_b = [d for d in a if d in b and d > 0]  # intersection of a and b
        assert not a_b, 'The two inputs should not have common elements'

        connected_dict = []

        # Find now which labels should be merged
        # mapped will be a list of 2d tuples. For example:
        # If mapped = [(7,2), (7,5), (8,1)]
        # it means that:
        #   label 7 and 2 should be considered the same
        #   label 7 and 5 should be considered the same, hence 2 and 5 are also the same
        #   label 8 and 1 should be considered the same
        t = set([d for d in list(zip(a, b)) if 0 not in d])
        mapped = list(zip(*t))

        if mapped:
            nlabels = np.array([a, b]).max()
            mat = _to_csr_matrix(mapped[0], mapped[1], nlabels + 1)
            n_components, labels = connected_components(csgraph=mat, directed=False, return_labels=True)
            connected_labels = _get_connected_labels(labels)

            _aa = []
            _bb = []
            for _list in connected_labels:
                _a = [lookup_label_a[d] for d in _list if d in a]
                _b = [lookup_label_b[d] for d in _list if d in b]

                connected_dict.append({'a': _a, 'b': _b})
        else:
            connected_labels = []
            connected_dict = []

        # print(connected_labels)
        return connected_dict

    def _shift_labels(self, a, b):
        # New labels must be negative
        #
        # Shifts the non-zero elements of a and b so that both lists a and b have values >= 0
        # then shifts a only if a and b have common elements so that they do not intersect
        # example: _shift_labels([2,1], [20,10]) gives:
        #           (array([2, 1]), array([20, 10]), {2: 2, 1: 1}, {20: 20, 10: 10})
        # nothing really changes since the input data are well formed, hence no shift of any kind has to be taken
        # For _shift_labels([-2,-1], [2,1]) then the output is
        #           (array([1, 2]), array([5, 4]), {1: -2, 2: -1}, {5: 2, 4: 1})
        # because [-2, -1] has be be shifted by 3 to become positive: [1, 2]. The same shift is also applied
        # to the second list, [2, 1], which becomes [5, 4]
        #
        a = np.array(a)
        b = np.array(b)

        _a = a.copy()
        _b = b.copy()

        mina = min(a)
        minb = min(b)

        if mina < 0 or minb < 0:
            k1 = abs(min(mina, minb)) + 1
        else:
            k1 = 0

        a[a != 0] = a[a != 0] + k1
        b[b != 0] = b[b != 0] + k1

        # take the intersection
        if np.any(np.in1d(a, b)):
            a_b = a[np.in1d(a, b)]
            a_b = a_b[a_b > 0]
        else:
            a_b = []

        if np.any(a_b) & (np.any(a) & np.any(b)):
            k2 = max([max(a), max(b)])
        else:
            k2 = 0

        a[a > 0] = a[a > 0] + k2

        # make a map to link the shifted labels with the original ones
        assert len(a) == len(_a)
        assert len(b) == len(_b)
        rmap_a = {a[i]: _a[i] for i, d in enumerate(a)}
        rmap_b = {b[i]: _b[i] for i, d in enumerate(b)}

        return a, b, rmap_a, rmap_b

    def collate_arrays(self, d):
        arr = self.tile_topo(d)
        stacked_rows = []
        for row in arr:
            row_temp = []
            for id in row:
                if np.isnan(id):
                    arr = np.zeros(self.tile_shape).astype(np.int32)
                else:
                    id = id.astype(np.int32)
                    arr = self.tiles[id]['label_image'].toarray().astype(np.int32)
                row_temp.append(arr)
            stacked_rows.append(np.hstack(row_temp))

        if len(stacked_rows) > 0:
            rows = self._padded(stacked_rows)  # <---- I THINK THIS IS NOT NEEDED ANYMORE
            return np.vstack(rows)
        else:
            return np.array([])

    def tile_topo(self, d):
        a = np.arange(self.tiles_down * self.tiles_across).reshape((self.tiles_down, self.tiles_across))
        mask = np.isin(a, d)
        return np.where(mask, a, np.nan)[mask.any(axis=1)][:, mask.any(axis=0)]

    def _padded(self, data):
        dims = max([d.shape for d in data])
        out = []
        for d in data:
            if d.shape != dims:
                p = np.zeros(dims)
                p[:d.shape[0], :d.shape[1]] = d
                out.append(p)
            else:
                out.append(d)
        return out

    def find_offset(self, tile_ids):
        sanity_check = np.array([self.tiles[d]['tile_id'] == d for d in tile_ids])
        assert np.all(sanity_check)
        offset_x = min([self.tiles[d]['tile_offset_x'] for d in tile_ids])
        offset_y = min([self.tiles[d]['tile_offset_y'] for d in tile_ids])
        return offset_x, offset_y

    def assign_cell_id(self):
        """ Add an extra column to be used as cell id
        This should be made redundant. The label can be used instead.
        """
        cell_id = self.cell_props.label - 1
        cell_id[cell_id < 0] = np.nan
        return cell_id

    def writer(self, dirpath):
        '''
        save the data to the flatfile
        :return:
        '''

        if os.path.exists(dirpath) and os.path.isdir(dirpath):
            shutil.rmtree(dirpath)
        os.mkdir(dirpath)

        # 1. save the cell props
        cell_props = self.cell_props.copy()
        cell_props['x'] = cell_props.x.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        cell_props['y'] = cell_props.y.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        cell_props['cell_id'] = cell_props.cell_id.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        cell_props['label'] = cell_props.label.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

        cells_headers = ['cell_id', 'label', 'tile_id', 'area', 'x', 'y']
        cell_props[cells_headers].to_csv(os.path.join(dirpath, '_cells.csv'), index=False)

        # 2. save the cell coords
        coords_headers = ['cell_id', 'label', 'coords']
        cell_props[coords_headers].to_json(os.path.join(dirpath, '_cellCoords.json'), orient='records')

        # 3. save the spots
        spots_df = self.spots.copy()
        spots_df['target'] = spots_df.Gene
        spots_df['x_global'] = spots_df.x
        spots_df['y_global'] = spots_df.y
        spots_df['tile_id'] = spots_df.tile_id
        spots_df['x_cell'] = spots_df.x_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        spots_df['y_cell'] = spots_df.y_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

        spots_headers = ['x_global', 'y_global', 'tile_id', 'label', 'target', 'x_cell', 'y_cell']
        spots_df[spots_headers].to_csv(os.path.join(dirpath, '_spots.csv'), index=False)
        logger.info('Total number of collected spots: %d' % spots_df.shape[0])

        return cell_props[cells_headers], cell_props[coords_headers], spots_df[spots_headers]

