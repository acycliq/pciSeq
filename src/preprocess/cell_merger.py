import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
import skimage.measure as skmeas
import itertools
from collections import defaultdict
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
from src.preprocess.cell_borders import extract_borders_par, get_label_contours
from src.preprocess.post import Post_merge
from src.preprocess.utils import _to_csr_matrix, _get_connected_labels
import logging

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



class Stage(object):
    def __init__(self, tile_obj, spots_all):
        self.counter = itertools.count()
        self.merge_register = Merge_register(self)
        self.tiles = tile_obj.tiles
        self.tiles_across = tile_obj.tiles_across
        self.tiles_down = tile_obj.tiles_down
        self.tile_shape = tile_obj.tile_shape
        self.scaling_factor = 1  # ti
        self.compose_dict(spots_all)
        # self._global_labels = None

    # @property
    # def global_labels(self):
    #     return self._global_labels

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
        """  Merge cells clipped by two or more fovs. """
        for tile in self.tiles:
            logger.info('\n')
            logger.info('Doing fov %i' % tile['tile_id'])
            self.merge(tile)
        logger.info('Relabelling finished')
        logger.info('\n')

    def merge(self, tile):
        """
        Does most of the heavy lifting for cell merging. Mutates in-place the label_image arrays of three fovs.
        If fov has fov_id = i then the mutated label_images are for the fovs with:
            fov_id = i
            fov_id = i + 1 (the neighbouring fov at the right)
            fov_id = i - #fovs_across (the neighbouring fov at the top)
        Parameters
        ----------
        fov: an instance of the class Fov

        Notes
        -----
        Is is assumed that each fov is big enough (relative to the cells) so that there is no cell bigger in size that a fov.
        For example, the most complicated case will be a cell clipped by four fovs forming a 2x2 setup with the cell centroid close
        at the intersection of the four fovs
        """
        tile_id = tile['tile_id']
        adj_img = self.adjacent_tile(tile_id)

        logger.info('fov_%d neighbours: (above, left): (%s, %s)' % (tile_id, adj_img['up'], adj_img['left']))

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
        Compares the label_image arrays from two neighbouring (one next another) fovs. If the last column of the
        label_image at the left and the first column of the one at the right have non-zero values at the same location
        then the labels at these locations are assigned a new and common label
        Parameters
        ----------
        adjc_fov: an instance of the class Fov
            The neighbouring fov. Could be the neighbour from the right, or from above
        fov: an instance of the class Fov
            the current fov
        transpose: bool. Optional
            if adjc_fov is the neighbour from the top, then set this to True. Default is False


        Returns
        -------
        temp_a, temp_b: tuple
            A tuple of two label_image arrays that correspond to the adjacent and the current fov respectively
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
                # logger.info('fov_%d: label %d ---> label %d' % (adjc_fov['fov_id'], x, new_label))

            for x in d['b']:
                temp_b.data[temp_b.data == x] = new_label
                self.merge_register.update_register(tile['tile_id'], new_label, x)
                # logger.info('fov_%d: label %d ---> label %d' % (fov['fov_id'], x, new_label))
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

    def calc_props(self):
        """ For the clipped cells, collate together all the necessary label_image arrays so that the cell is
        not clipped anymore. Then read that mosaic of arrays and get cell centroid, area and other properties.
        For the unclipped cells these properties are kept in image_objects.
        Combine all (clipped and unclipped) and make a single dataframe with columns: 'label', 'fov_id', 'area',
        'x_local', 'y_local', 'fov_offset_x', 'fov_offset_y', 'is_clipped', 'x', 'y', 'coords'.
        It also assigns a unique global label to the cells. The last column with label 'coords' keeps the
        coordinates of the cell boundary
        """

        # make a dict where each key is a list of fovs and each value is a list of labels. The latter
        # holds the labels whose corresponding objects as these were identified by the cell segmentation
        # extend outside one single fov and span across the fovs kept in the corresponding key of the dict
        v = defaultdict(list)
        for key, value in sorted(list(self.merge_register.entries.items())):
            if key not in v[tuple(value)]:
                v[tuple(value)].append(key)

        # first, do the merge cells
        # check the register to see which cells are merged
        merged_props_dict = {'label': [], 'tile_id': [], 'area': [], 'x_local': [], 'y_local': []}
        for tt in v.items():
            labels = tt[1]
            tile_ids = sorted(tt[0])

            logger.info('calculating centroid and area for the merged cells with labels %s' % labels)
            if len(tile_ids) > 1:
                label_image = self.collate_arrays(tile_ids)
            else:
                label_image = None  # you shouldnt get here. Since we merge cells we need at least two fovs

            # sanity check
            coo_mat = set(coo_matrix(label_image).data)
            assert np.all([label in coo_mat for label in labels]), 'a label that should exist in the collated label_image seem to be missing'

            props = skmeas.regionprops(label_image.astype(np.int32))
            clipped_cells = list(filter(lambda d: d.label in labels, props))
            # unclipped_cells = list(filter(lambda d: d.label not in labels, props))
            for i, p in enumerate(clipped_cells):
                centroid_tile_id, centroid_coords = self.locate_tile(p.centroid, tile_ids)

                logger.info('cell with label %d is clipped by fov_ids %s' % (p.label, tile_ids))

                # # sanity check
                # filtered_fovs = list(filter(lambda d: d['fov_id'] in fov_ids, self.fovs))
                # mask = [d['image_objects'].label == labels[i] for d in filtered_fovs]
                # total_area = sum([d['image_objects'].area[mask[i]].values[0] for i, d in enumerate(filtered_fovs)])
                # # logger.info('total sum for the sanity check is: %d' % total_area)
                # assert int(total_area) == int(p.area)

                merged_props_dict['label'].append(p.label)
                merged_props_dict['tile_id'].append(centroid_tile_id)
                merged_props_dict['area'].append(p.area)
                merged_props_dict['x_local'].append(centroid_coords[1])
                merged_props_dict['y_local'].append(centroid_coords[0])

                # sort the dict by the label
                merged_props_dict = pd.DataFrame(merged_props_dict).sort_values(by=['label']).to_dict(orient='list')
                # logger.info('merged_props_dict ends')


        img_obj = pd.concat([d['image_objects'] for d in self.tiles])
        dup_vals = list(set(img_obj[img_obj.duplicated('label')].label.values))
        # np.setdiff1d(np.array(self.cell_props['label']), dup_vals) # <--  WHY I HAVE THAT?

        unclipped_labels = list(set(img_obj.label.values) - set(merged_props_dict['label']))
        unclipped_cells = img_obj.iloc[np.isin(img_obj.label.values, unclipped_labels), :]

        # sanity check. Cannot have duplicate labels. If a label is duplicate then the cell is clipped, hence it should have been captured by
        # 'merged_props_dict'
        assert unclipped_cells[unclipped_cells.duplicated('label')].empty, "Dataframe 'unclipped_img_obj' contains dublicate labels."

        # unclipped_cells = img_obj.drop_duplicates('label', keep=False)
        merged = pd.DataFrame(merged_props_dict)

        df_1 = unclipped_cells.merge(pd.DataFrame(self.tiles)[['tile_id', 'tile_offset_x', 'tile_offset_y']], how='left', on='tile_id')
        df_1['is_clipped'] = False
        df_2 = merged.merge(pd.DataFrame(self.tiles)[['tile_id', 'tile_offset_x', 'tile_offset_y']], how='left', on='tile_id')
        df_2['is_clipped'] = True

        cell_props = pd.concat([df_1.reset_index(drop=True), df_2.reset_index(drop=True)], axis=0).reset_index(drop=True)
        cell_props['x'] = cell_props.x_local + cell_props.tile_offset_x
        cell_props['y'] = cell_props.y_local + cell_props.tile_offset_y

        # sort by x and y and then relabel so that labelling looks more deterministic.
        # Not necessary, it is nice to have tho. If however that takes place then the
        # label_image array has to be updated to be inline with the new labels
        cell_props = cell_props.sort_values(by=['x', 'y'])

        global_label = np.arange(cell_props.shape[0]) + 1
        label_map = np.column_stack((cell_props.label.values, global_label))
        cell_props['label'] = global_label
        cell_props['label'] = global_label

        logger.info('')
        logger.info('updating the labels on label_image to align them with cell_props')
        for tile in self.tiles:
            data = self.update_label_image(tile, label_map)
            tile['label_image'].data = data.astype(np.int64)

        logger.info('Done')

        # relabel now the key in the merge_register dict to align them
        # with the global labels
        self.merge_register = self._update_dict_keys(label_map)

        # DO I ALSO NEED TO UPDATE THE IMAGE_OBJECTS FOR EACH FOV??
        # Labels there need to be updated for consistency
        # However in a cell is split in to two parts and exists in two
        # fovs then x_local and y_local coords show the centroid of the partial shape
        # The array cell_props is the correct place for such lookups.

        # Find now the outline of the cells
        _cell_boundaries = self.cell_boundaries(cell_props)
        out = cell_props.merge(_cell_boundaries, how='left', on=['label'])

        return out


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

    def locate_tile(self, centroid, tile_ids):
        # first get the fov_id
        t = self.tile_topo(tile_ids)
        row = int(centroid[0] // self.tile_shape[0])  # <-- I think I should be dividing by fov_shape[1] instead
        col = int(centroid[1] // self.tile_shape[1])  # <-- I think I should be dividing by fov_shape[0] instead

        tile_id = t[row, col]
        assert ~np.isnan(tile_id)
        assert tile_id in tile_ids

        # calc the local coordinates
        coord_row = centroid[0] % self.tile_shape[0]  # <-- I think I should be dividing by fov_shape[1] instead
        coord_col = centroid[1] % self.tile_shape[1]  # <-- I think I should be dividing by fov_shape[0] instead

        return tile_id, (coord_row, coord_col)

    def update_label_image(self, tile, label_map):
        _x = tile['label_image'].data
        _y = label_map[:, 0]
        idx = self._remap(_x, _y)
        return label_map[idx, 1]

    def _remap(self, x, y):
        index = np.argsort(y)
        sorted_y = y[index]
        sorted_index = np.searchsorted(sorted_y, x)

        xindex = np.take(index, sorted_index, mode="clip")
        mask = y[xindex] != x

        result = np.ma.array(xindex, mask=mask)
        return result

    def _update_dict_keys(self, label_map):
        # relabel now the key in the merge_register dict to align them
        # with the global labels
        dict_keys = np.array(list(self.merge_register.entries.keys()))
        keys_mask = self._remap(dict_keys, label_map[:, 0])
        global_keys = label_map[keys_mask, 1]
        dict_vals = list(self.merge_register.entries.values())
        global_dict = dict(zip(global_keys, dict_vals))
        return global_dict



    def cell_boundaries(self, cell_props):
        '''
        calculate the outlines of the cells
        :return:
        '''

        # loop over the self.cell_props
        res_list = []
        for tile in self.tiles:
            if np.any(tile['label_image'].data):
                df = self.obj_outline(tile, cell_props)
                res_list.append(df)
            else:
                logger.info('tile:%d empty, No cells to draw boundaries were found' % tile['tile_id'])
        _df = pd.concat(res_list).astype({"label": int})

        # make a Dataframe to keep boundaries of the cells which are not clipped by the fov
        df_1 = _df.iloc[np.isin(_df.label, cell_props[~cell_props.is_clipped].label)]

        # get the labels of the clipped cells
        in_multiple_tiles = sorted(cell_props[cell_props.is_clipped].label.values)
        logger.info('There are %d cells whose boundaries span across multiple fovs' % len(in_multiple_tiles))

        # find the boundaries of the clipped cells
        _list = self.collate_borders_par(in_multiple_tiles)
        df_2 = pd.DataFrame(_list).astype({"label": int})

        # Both clipped and unclipped in a dataframe
        res = pd.concat([df_1, df_2])

        set_diff = set(cell_props.label) - set(res.label.values)
        if set_diff:
            unresolved_labels = pd.DataFrame({'label': list(set_diff), 'coords': np.nan * np.ones(len(set_diff))})
            res = pd.concat([res, unresolved_labels])

        # assert set(_df.label.values) == set(res.label.values)
        assert res.shape[0] == cell_props.shape[0]
        assert np.all(sorted(res.label) == sorted(cell_props.label))
        assert np.unique(res.label).size == res.shape[0], 'Array cannot have duplicates'
        return res.sort_values(['label'], ascending=[True])


    def collate_borders_par(self, in_multiple_fovs):
        n = max(1, cpu_count() - 1)
        pool = ThreadPool(16)
        results = pool.map(self.collate_borders_helper, in_multiple_fovs)
        pool.close()
        pool.join()
        return results

    def collate_borders_helper(self, label):
        out = {}
        logger.info('label: %d. Finding the cell boundaries' % label)
        label_image = self.collate_arrays(self.merge_register[label])
        offset_x, offset_y = self.find_offset(self.merge_register[label])
        out['coords'] = get_label_contours(label_image, label, offset_x, offset_y)
        out['label'] = label

        # ALTERNATIVE:
        # Maybe this is better since it using the same code as the function
        # that calcs non-clipped cells.
        # Needs futher testing, but the code will be something close to this:
        #
        # logger.info('Using chain codes')
        # exclude = set(label_image[label_image != label])
        # temp = extract_borders_dip(label_image, offset_x, offset_y, exclude)
        # logger.info('Cell boundaries (chain codes) found')
        return out

    def find_offset(self, tile_ids):
        sanity_check = np.array([self.tiles[d]['tile_id'] == d for d in tile_ids])
        assert np.all(sanity_check)
        offset_x = min([self.tiles[d]['tile_offset_x'] for d in tile_ids])
        offset_y = min([self.tiles[d]['tile_offset_y'] for d in tile_ids])
        return offset_x, offset_y

    def obj_outline(self, tile, cell_props):
        logger.info('Getting cell boundaries for cells in tile: %d' % tile['tile_id'])
        label_image = tile['label_image'].toarray()
        offset_x = tile['tile_offset_x']
        offset_y = tile['tile_offset_y']
        clipped_cells = cell_props[cell_props.is_clipped].label.values

        df = extract_borders_par(label_image, offset_x, offset_y, clipped_cells)
        temp_df = extract_borders_par(label_image, offset_x, offset_y, clipped_cells)
        return df

    def assign_cell_id(self):
        """ Add an extra column to be used as cell id
        This should be made redundant. The label can be used instead.
        """
        cell_id = self.cell_props.label - 1
        cell_id[cell_id < 0] = np.nan
        return cell_id

    def writer(self):
        '''
        save the data to the flatfile
        :return:
        '''

        # 1. save the cell props
        cell_props = self.cell_props.copy()
        cell_props['x'] = cell_props.x.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        cell_props['y'] = cell_props.y.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        cell_props['cell_id'] = cell_props.cell_id.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        cell_props['label'] = cell_props.label.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

        cells_headers = ['cell_id', 'label', 'tile_id', 'area', 'x', 'y']
        cell_props[cells_headers].to_csv('cells.csv', index=False)

        # 2. save the cell coords
        coords_headers = ['cell_id', 'label', 'coords']
        cell_props[coords_headers].to_json('cell_coords.json', orient='records')

        # 3. save the spots
        spots_df = self.spots.copy()
        spots_df['target'] = spots_df.Gene
        spots_df['x_global'] = spots_df.x
        spots_df['y_global'] = spots_df.y
        spots_df['fov_id'] = spots_df.fov_id
        spots_df['x_cell'] = spots_df.x_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        spots_df['y_cell'] = spots_df.y_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

        spots_headers = ['x_global', 'y_global', 'fov_id', 'label', 'target', 'x_cell', 'y_cell']
        spots_df[spots_headers].to_csv('spots.csv', index=False)
        logger.info('Total number of collected spots: %d' % spots_df.shape[0])

