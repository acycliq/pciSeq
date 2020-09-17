import numpy as np
import pandas as pd
import urllib, base64
import json
import os
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
from scipy import ndimage
import config
import skimage.measure as skmeas
import itertools
# import credentials
import math
from collections import defaultdict
from src.preprocess.fov import Fov
# from preprocess import Stage_spots
from multiprocessing import Pool as ProcessPool
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
from src.preprocess.cell_borders import extract_borders_par, get_label_contours
from src.preprocess.cell_borders import extract_borders
import logging
import time

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)

logger = logging.getLogger()


def _to_csr_matrix(i, j, n):
    """Using i and j as coo-format coordinates, return csr matrix."""
    n = int(n)
    v = np.ones_like(i)
    mat = coo_matrix((v, (i, j)), shape=(n, n))
    return mat.tocsr()


def get_dir(my_config, fov_id):
    root_dir = my_config['FOV_ROOT']
    return os.path.join(root_dir, 'fov_' + str(fov_id))


def _get_connected_labels(mylist):
    '''
    find which positions in the input list have repeated values
    Example:
     If mylist = [0, 4, 4] then it returns [1, 2] because 4 appears twice, at positions 1 and 2 of the input list
     if mylist = [0,1,2,1,3,4,2,2] then it returns [[1, 3], [2, 6, 7]]
    :param mylist:
    :return:
    '''
    output = defaultdict(list)
    # Loop once over mylist, store the indices of all unique elements
    for i, el in enumerate(mylist):
        output[el].append(i)
    return np.array([np.array(d) for d in output.values() if len(d) > 1])


class Stage(object):
    def __init__(self, fovs_obj, spots_all, cellmaps):
        """
        The Stage object contains the tools to tell if a spot falls within the cell soma (and which cell that is),
        merge cells that span over many fovs and also label these cells with a unique id

        Parameters
        ----------
        fovs_obj: an object of type Fov
        spots_all: a dataframe with columns ['Gene', 'x', 'y']. It contains all the spots. Column ```Gene``` is the gene name
                    and columns '''x''', '''y''' contain the x,y coordinated of the cell centroid respectively.
        cellmaps: a list containing the ```label_image''' arrays for each fov. A '''label_image''' is an array where the non-zero
                    values denote the label of the physical object (ie cell) that correspond to that specific location. Zero means
                    that the corresponding pixel falls on the background

        Attributes
        ----------
        my_counter: int
            counter that is called internally while new labels to the cells
        merge_register: dict
            dictionary that keeps the newly created labels (as keys) and the fovs where these labels appear (as values)
        cell_props: dataframe
            keeps some useful cell properties
        fovs: list
            Each element of the list is an oject of type FOV
        fovs_across: int
            number of fovs along the x-axis
        fovs_down: int
            number of fovs along the y-axis
        fov_shape: list
            list with 2 elements denoting the length of the x-side and y-side of the fov. All fovs are assumed to be of the
            same shape


        Note: '''label_image''' is typically obtained via cell segmentation. The '''cellmaps''' list contains the output
                after cell-segmenting each fov that composes our big dapi image we performed the experiment on

        """
        self.my_counter = itertools.count()
        self.merge_register = defaultdict(list)
        self.cell_props = None
        self.spots = None
        self.cellmaps = cellmaps
        self.fovs = fovs_obj.fovs
        self.fovs_across = fovs_obj.fovs_across
        self.fovs_down = fovs_obj.fovs_down
        self.fov_shape = fovs_obj.fov_shape
        self.scaling_factor = fovs_obj.scaling_factor
        self.compose_dict(spots_all)

    def compose_dict(self, spots_all):
        """
        Mutates in-place the fov object by adding two more key/value pairs

        Parameters
        ----------
        spots_all: dataframe
            Contains all the spots (Gene names and x,y coords) for the full image
        """
        for i, d in enumerate(self.fovs):
            d['label_image'] = self.fov_label_image(i)  # label_image for the i-th fov
            d['spots'] = self.fov_spots(spots_all, i)   # spots for the i-th fov

    def fov_label_image(self, i):
        """ label_image for the i-th fov """
        if self.cellmaps is not None:
            label_image = coo_matrix(self.cellmaps[i])
        else:
            raise ValueError('cellmaps is not defined')
        return label_image

    def fov_spots(self, data, i):
        """ spots for the i-th fov """
        x_range = self.fovs[i]['fov_range']['x']
        y_range = self.fovs[i]['fov_range']['y']
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

    def assign_spot_parent(self):
        res_list = []
        for i, fov in enumerate(self.fovs):
            assert i == fov['fov_id'], 'The list is miss-aligned'
            # Assign the parent cell and its coordinates
            res = self._assign_spot_parent_helper(fov)
            fov['spots'] = res
            res['fov_id'] = fov['fov_id']
            res_list.append(res)
        return pd.concat(res_list).astype({"label": int})

    def _assign_spot_parent_helper(self, fov):
        '''
        assign the parent cell
        :return:
        '''

        # first find the label of the parent cell
        spots_df = fov['spots']
        label_image = fov['label_image']
        spots_temp = self.spot_parent_label(fov, spots_df.copy(), label_image)

        # find now the coordinates of the parent cell
        out = self.spot_parent_coords(spots_temp.copy())
        return out

    def spot_parent_label(self, fov, spots, sp_label_image):
        # spots = fov['spots']
        if sp_label_image.nnz == 0:
            # The tile empty, all spots are on the background
            spots['label'] = 0
        else:
            # 1. get the label_image
            label_image = sp_label_image.toarray()

            # 2. unscaled and local coordinates for the spots (as a sparse array)
            coo = self.coofy(spots[['x', 'y']].copy(), fov)

            coo_arr = coo.toarray()
            label_coords = coo_matrix(label_image * coo_arr.astype(bool))
            spot_id_coords = coo_matrix(label_image.astype(bool) * coo_arr)

            df = pd.DataFrame({'x': label_coords.col,
                               'y': label_coords.row,
                               'label': label_coords.data},
                              index=spot_id_coords.data).astype(int)

            if np.any(df.index.duplicated()):
                logger.warning('Found %d duplicated. Investigate!' % df.index.duplicated().sum())
                ## that actually means that the same spot (ie read/dot) exists at two different locations at the same time

            df = df[~df.index.duplicated()]
            spots['label'] = df.label

            # if nan it means the spot is on the background. Hence set the label = 0
            spots['label'] = spots.label.fillna(0).astype(int)
        return spots

    def spot_parent_coords(self, spots):
        # spots = fov['spots']
        cells = self.cell_props[['label', 'x', 'y']].rename(columns={'x': 'x_cell', 'y': 'y_cell'})
        out = spots.merge(cells, how='left', on=['label'])
        out['label'] = out.label.astype('Int64')
        return out

    def localise_coords(self, spots, fov):
        '''
        convert spots coords to local coords
        (lacal means the origin is the top left corner of the fov not of the full image)
        :param spots:
        :return:
        '''
        # spots = fov['spots']
        df = spots[['x', 'y']]
        x0 = fov['fov_offset_x']
        y0 = fov['fov_offset_y']
        origin = np.array([x0, y0])
        df = df - origin
        return df.astype(int)

    def coofy(self, spots, fov):
        spots_loc = self.localise_coords(spots, fov)
        x = spots_loc.x.values.astype(int)
        y = spots_loc.y.values.astype(int)

        # make a sparse array of size 6000-6000
        # _len = np.diff(np.array(fov['fov_range']['x'])) * fov['fov_scaling_factor']
        # _height = np.diff(np.array(fov['fov_range']['y'])) * fov['fov_scaling_factor']

        # coo_shape = self.fov_shape + 1
        coo_shape = np.array(fov['label_image'].shape) + 1
        coo = coo_matrix((spots.index.values, (y, x)), shape=coo_shape)
        # Why the +1? starfish includes the boundaries,
        # hence the extra pixel. That sounds like a bug maybe? a tile of size 2000-by-2000 should have a range
        # [0, 1999] going across and [0, 1999] going down

        # ok, that a bit of a mystery!
        coo_arr = coo.toarray()
        out = coo_matrix(coo_arr[:-1, :-1])  # why the -1? See comment right above
        return out

    def adjacent_fov(self, fov_id):
        if fov_id % self.fovs_across != 0:
            left = fov_id - 1
        else:
            left = None

        if fov_id >= self.fovs_across:
            up = fov_id - self.fovs_across
        else:
            up = None
        return {'left': left, 'up': up}

    def merge_cells(self):
        """  Merge cells clipped by two or more fovs. """
        for fov in self.fovs:
            logger.info('\n')
            logger.info('Doing fov %i' % fov['fov_id'])
            self.merge(fov)
        logger.info('Relabelling finished')
        logger.info('\n')

    def merge(self, fov):
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
        fov_id = fov['fov_id']
        adj_img = self.adjacent_fov(fov_id)

        logger.info('fov_%d neighbours: (above, left): (%s, %s)' % (fov_id, adj_img['up'], adj_img['left']))

        # keep the original label_image somewhere
        fov['label_image_zero'] = fov['label_image']

        # Bottom border of the label array above
        if (adj_img['up'] is not None) and np.any(self.fovs[adj_img['up']]['label_image'].data):
            fov_up = self.fovs[adj_img['up']]
            coo_aa, coo_bb = self.dissolve_borders(fov_up, fov, transpose=True)
            fov_up['label_image'] = coo_aa
            fov['label_image'] = coo_bb

        if adj_img['left'] is not None:
            fov_left = self.fovs[adj_img['left']]
            coo_a, coo_b = self.dissolve_borders(fov_left, fov)
            fov_left['label_image'] = coo_a
            fov['label_image'] = coo_b

    def dissolve_borders(self, adjc_fov, fov, transpose=False):
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
            adjc_img = adjc_fov['label_image'].transpose()
            img = fov['label_image'].transpose()
        else:
            adjc_img = adjc_fov['label_image']
            img = fov['label_image']

        arr = adjc_img.toarray()
        adjc_border = arr[:, -1]
        border = img.toarray()[:, 0]

        logger.info('length of adjc_border: %d' % adjc_border.shape[0])
        logger.info('length of adjc_border: %d' % border.shape[0])
        matched_labels = self.connect_labels(adjc_border, border)
        temp_a = self.fovs[adjc_fov['fov_id']]['label_image'].copy()
        temp_b = self.fovs[fov['fov_id']]['label_image'].copy()

        for d in matched_labels:
            new_label = self._new_label(d)
            for x in d['a']:
                temp_a.data[temp_a.data == x] = new_label
                self.update_register(adjc_fov['fov_id'], new_label, x)
                # logger.info('fov_%d: label %d ---> label %d' % (adjc_fov['fov_id'], x, new_label))

            for x in d['b']:
                temp_b.data[temp_b.data == x] = new_label
                self.update_register(fov['fov_id'], new_label, x)
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

    def label_generator(self):
        return -1 * (next(self.my_counter) + 1)

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

    def image_objects(self, fov):
        """
        Reads a label_image and returns a dataframe with the objects found on it.
        The array has index the label and columns: fov_id, area and the cell centroids x, y
        """
        fov_id = fov['fov_id']
        label_image = fov['label_image'].toarray().astype(np.int32)

        props = skmeas.regionprops(label_image)
        props_arr = [[int(p.label), int(fov_id), int(p.area), p.centroid[1], p.centroid[0]] for i, p in
                     enumerate(props)]
        props_df = pd.DataFrame(props_arr, columns=['label', 'fov_id', 'area', 'x_local', 'y_local'])

        # stick the df to the dict
        fov['image_objects'] = props_df  # <--- Not sure If I should do this or just return the df and add the dictionary filed at the calliing function
        return fov

    def global_labels(self):
        """ Assign new labels to the cells. I will also create a new key '''image_objects''' in the fov dict
        with value a dataframe that keeps the label, fov_id, area, and centroid (local) coords of the cells
        """
        self.global_labels_par()
        # Need also to flip the sign of the labels in the merge_register dict
        self.flip_sign()

    def global_labels_par(self):
        n = max(1, cpu_count() - 1)
        pool = ThreadPool(n)
        results = pool.map(self.global_labels_helper, self.fovs)
        pool.close()
        pool.join()

        return results

    # def global_labels(self):
    #     '''
    #     DEPRECATED: Replaced the function that uses parallelism
    #     Same as "global_labels_par()" but without parallelism
    #     :return:
    #     '''
    #     for fov in self.fovs:
    #         # logger.info('fov:%d, Setting global labels' % fov['fov_id'])
    #         data = fov['label_image'].data
    #         labels = np.unique(data[data > 0])
    #         label_map = {d: self.label_generator() for d in labels}
    #         fov['label_image'].data = np.array([label_map[d] if d > 0 else d for d in data])
    #         logger.info('fov: %d: label map is: %s' % (fov['fov_id'], label_map))
    #
    #         clipped_cell_labels = {k: v for k, v in self.merge_register.items() if fov['fov_id'] in v}.keys()
    #         assert np.all([d in fov['label_image'].data for d in clipped_cell_labels]), \
    #             "A label that the cell register says should exist in the fov %d doesnt seem to be verified by the label_image" % \
    #             fov['fov_id']
    #
    #         fov['label_image'].data = -1 * fov['label_image'].data
    #         if np.any(fov['label_image'].data):
    #             logger.info('fov:%d, Global labels are set.' % fov['fov_id'])
    #         else:
    #             logger.info('fov:%d empty, nothing to do here' % fov['fov_id'])
    #
    #         self.image_objects(fov)
    #
    #     # finally flip the sign on the dict
    #     self.flip_sign()


    def global_labels_helper(self, fov):
        logger.info('fod_id is: %d' % fov['fov_id'])
        data = fov['label_image'].data

        # Get the positive labels. Negative means the cell is clipped, hence it has been
        # given a new label (negative valued)
        labels = np.unique(data[data > 0])
        label_map = {d: self.label_generator() for d in labels}
        fov['label_image'].data = np.array([label_map[d] if d > 0 else d for d in data])

        # Sanity check. merge register keeps the cells that are clipped, the cells that span
        # over more than a single fov
        # Get the labels from the merge register for this particular fov
        clipped_cell_labels = {k: v for k, v in self.merge_register.items() if fov['fov_id'] in v}.keys()
        assert np.all([d in fov['label_image'].data for d in clipped_cell_labels]), \
            "A label that the cell register says should exist in the fov %d doesnt seem to be verified by the " \
            "label_image" % fov['fov_id']

        # Flip now the sign of the labels
        fov['label_image'].data = -1 * fov['label_image'].data
        if np.any(fov['label_image'].data):
            logger.info('fov:%d, Global labels are set.' % fov['fov_id'])
        else:
            logger.info('fov:%d empty, nothing to do here' % fov['fov_id'])

        # Finally, get area, centroid coords and label of the cells you find on the fov
        self.image_objects(fov)
        return True

    def flip_sign(self):
        # need to flip the sign of the keys in the merge_register dict
        # Maybe I should add some sanity checks here.
        # For example I can do the flipping inside the loop of the calling function and while
        # you loop get the keys (ie the labels) for the particular fov and then flip the label
        # sign (if it negative because it may have been flipped already if the label exists in
        # more that one fov

        if bool(self.merge_register):
            my_keys = list(self.merge_register.keys())
            vals = list(self.merge_register.values())

            assert max(my_keys) < 0, 'All keys should be non-zero negative'
            new_keys = [-1 * d for d in my_keys]

            self.merge_register = dict(zip(new_keys, vals))
        else:
            logger.info('Merge register empty')

    def update_register(self, fov_id, label, old_label):
        self.merge_register[label].append(fov_id)
        self.merge_register[label] = sorted(list(set(self.merge_register[label])))

        logger.info('fov_%d: label %d ---> label %d' % (fov_id, old_label, label))
        if (old_label != label) and (old_label < 0):
            self.replace_label(old_label, label)

    def replace_label(self, old_label, label):
        # replace the register
        _dict = self.merge_register
        for fov in _dict[old_label]:
            if fov not in _dict[label]:
                _dict[label].append(fov)
            mask = self.fovs[fov]['label_image'].data == old_label
            self.fovs[fov]['label_image'].data[mask] = label
            logger.info('fov: %d: replaced labels in "label_image" that were equal to %d with %d' % (fov, old_label, label))

        logger.info('Dropped key, value pair: (%d, %s) from the merge_register' % (old_label, _dict[old_label]))
        _dict.pop(old_label)


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
        for key, value in sorted(list(self.merge_register.items())):
            if key not in v[tuple(value)]:
                v[tuple(value)].append(key)

        # first, do the merge cells
        # check the register to see which cells are merged
        merged_props_dict = {'label': [], 'fov_id': [], 'area': [], 'x_local': [], 'y_local': []}
        for tt in v.items():
            labels = tt[1]
            fov_ids = sorted(tt[0])

            logger.info('calculating centroid and area for the merged cells with labels %s' % labels)
            if len(fov_ids) > 1:
                label_image = self.collate_arrays(fov_ids)
            else:
                label_image = None  # you shouldnt get here. Since we merge cells we need at least two fovs

            # sanity check
            coo_mat = set(coo_matrix(label_image).data)
            assert np.all([label in coo_mat for label in labels]), 'a label that should exist in the collated label_image seem to be missing'

            props = skmeas.regionprops(label_image.astype(np.int32))
            clipped_cells = list(filter(lambda d: d.label in labels, props))
            # unclipped_cells = list(filter(lambda d: d.label not in labels, props))
            for i, p in enumerate(clipped_cells):
                centroid_fov_id, centroid_coords = self.locate_fov(p.centroid, fov_ids)

                logger.info('cell with label %d is clipped by fov_ids %s' % (p.label, fov_ids))

                # # sanity check
                # filtered_fovs = list(filter(lambda d: d['fov_id'] in fov_ids, self.fovs))
                # mask = [d['image_objects'].label == labels[i] for d in filtered_fovs]
                # total_area = sum([d['image_objects'].area[mask[i]].values[0] for i, d in enumerate(filtered_fovs)])
                # # logger.info('total sum for the sanity check is: %d' % total_area)
                # assert int(total_area) == int(p.area)

                merged_props_dict['label'].append(p.label)
                merged_props_dict['fov_id'].append(centroid_fov_id)
                merged_props_dict['area'].append(p.area)
                merged_props_dict['x_local'].append(centroid_coords[1])
                merged_props_dict['y_local'].append(centroid_coords[0])

                # sort the dict by the label
                merged_props_dict = pd.DataFrame(merged_props_dict).sort_values(by=['label']).to_dict(orient='list')
                # logger.info('merged_props_dict ends')


        img_obj = pd.concat([d['image_objects'] for d in self.fovs])
        dup_vals = list(set(img_obj[img_obj.duplicated('label')].label.values))
        # np.setdiff1d(np.array(self.cell_props['label']), dup_vals) # <--  WHY I HAVE THAT?

        unclipped_labels = list(set(img_obj.label.values) - set(merged_props_dict['label']))
        unclipped_cells = img_obj.iloc[np.isin(img_obj.label.values, unclipped_labels), :]

        # sanity check. Cannot have duplicate labels. If a label is duplicate then the cell is clipped, hence it should have been captured by
        # 'merged_props_dict'
        assert unclipped_cells[unclipped_cells.duplicated('label')].empty, "Dataframe 'unclipped_img_obj' contains dublicate labels."

        # unclipped_cells = img_obj.drop_duplicates('label', keep=False)
        merged = pd.DataFrame(merged_props_dict)

        df_1 = unclipped_cells.merge(pd.DataFrame(self.fovs)[['fov_id', 'fov_offset_x', 'fov_offset_y']], how='left', on='fov_id')
        df_1['is_clipped'] = False
        df_2 = merged.merge(pd.DataFrame(self.fovs)[['fov_id', 'fov_offset_x', 'fov_offset_y']], how='left', on='fov_id')
        df_2['is_clipped'] = True

        cell_props = pd.concat([df_1.reset_index(drop=True), df_2.reset_index(drop=True)], axis=0).reset_index(drop=True)
        cell_props['x'] = cell_props.x_local + cell_props.fov_offset_x
        cell_props['y'] = cell_props.y_local + cell_props.fov_offset_y

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
        for fov in self.fovs:
            data = self.update_label_image(fov, label_map)
            fov['label_image'].data = data.astype(np.int64)

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

    def _update_dict_keys(self, label_map):
        # relabel now the key in the merge_register dict to align them
        # with the global labels
        dict_keys = np.array(list(self.merge_register.keys()))
        keys_mask = self.my_remap(dict_keys, label_map[:, 0])
        global_keys = label_map[keys_mask, 1]
        dict_vals = list(self.merge_register.values())
        global_dict = dict(zip(global_keys, dict_vals))
        return global_dict

    def update_label_image(self, fov, label_map):
        _x = fov['label_image'].data
        _y = label_map[:, 0]
        idx = self.my_remap(_x, _y)
        return label_map[idx, 1]

    def my_remap(self, x, y):
        index = np.argsort(y)
        sorted_y = y[index]
        sorted_index = np.searchsorted(sorted_y, x)

        xindex = np.take(index, sorted_index, mode="clip")
        mask = y[xindex] != x

        result = np.ma.array(xindex, mask=mask)
        return result

    def locate_fov(self, centroid, fov_ids):
        # first get the fov_id
        t = self.fov_topo(fov_ids)
        row = int(centroid[0] // self.fov_shape[0])  # <-- I think I should be dividing by fov_shape[1] instead
        col = int(centroid[1] // self.fov_shape[1])  # <-- I think I should be dividing by fov_shape[0] instead

        fov_id = t[row, col]
        assert ~np.isnan(fov_id)
        assert fov_id in fov_ids

        # calc the local coordinates
        coord_row = centroid[0] % self.fov_shape[0]  # <-- I think I should be dividing by fov_shape[1] instead
        coord_col = centroid[1] % self.fov_shape[1]  # <-- I think I should be dividing by fov_shape[0] instead

        return fov_id, (coord_row, coord_col)

    def in_fov(self, centroid, fov):
        x_range = np.array(fov['fov_range']['x'])
        y_range = np.array(fov['fov_range']['y'])
        x = centroid[0]
        y = centroid[1]
        if (x_range[1] >= x >= x_range[0]) and (y_range[1] >= y >= y_range[0]):
            return True
        else:
            return False

    def collate_arrays(self, d):
        arr = self.fov_topo(d)
        stacked_rows = []
        for row in arr:
            row_temp = []
            for id in row:
                if np.isnan(id):
                    arr = np.zeros(self.fov_shape).astype(np.int32)
                else:
                    id = id.astype(np.int32)
                    arr = self.fovs[id]['label_image'].toarray().astype(np.int32)
                row_temp.append(arr)
            stacked_rows.append(np.hstack(row_temp))

        if len(stacked_rows) > 0:
            rows = self._padded(stacked_rows)  # <---- I THINK THIS IS NOT NEEDED ANYMORE
            return np.vstack(rows)
        else:
            return np.array([])

    def fov_topo(self, d):
        a = np.arange(self.fovs_down * self.fovs_across).reshape((self.fovs_down, self.fovs_across))
        mask = np.isin(a, d)
        return np.where(mask, a, np.nan)[mask.any(axis=1)][:, mask.any(axis=0)]

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

        cells_headers = ['cell_id', 'label', 'fov_id', 'area', 'x', 'y']
        cell_props[cells_headers].to_csv('expanded_cells_david_2.csv', index=False)

        # 2. save the cell coords
        coords_headers = ['cell_id', 'label', 'coords']
        cell_props[coords_headers].to_json('cell_coords_david_2.json', orient='records')

        # 3. save the spots
        spots_df = self.spots.copy()
        spots_df['target'] = spots_df.Gene
        spots_df['x_global'] = spots_df.x
        spots_df['y_global'] = spots_df.y
        spots_df['fov_id'] = spots_df.fov_id
        spots_df['x_cell'] = spots_df.x_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        spots_df['y_cell'] = spots_df.y_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

        spots_headers = ['x_global', 'y_global', 'fov_id', 'label', 'target', 'x_cell', 'y_cell']
        spots_df[spots_headers].to_csv('spots_david_2.csv', index=False)
        logger.info('Total number of collected spots: %d' % spots_df.shape[0])


    def cell_boundaries(self, cell_props):
        '''
        calculate the outlines of the cells
        :return:
        '''

        # loop over the self.cell_props
        res_list = []
        for fov in self.fovs:
            if np.any(fov['label_image'].data):
                df = self.obj_outline(fov, cell_props)
                res_list.append(df)
            else:
                logger.info('fov:%d empty, No cells to draw boundaries were found' % fov['fov_id'])
        _df = pd.concat(res_list).astype({"label": int})

        # make a Dataframe to keep boundaries of the cells which are not clipped by the fov
        df_1 = _df.iloc[np.isin(_df.label, cell_props[~cell_props.is_clipped].label)]

        # get the labels of the clipped cells
        in_multiple_fovs = sorted(cell_props[cell_props.is_clipped].label.values)
        logger.info('There are %d cells whose boundaries span across multiple fovs' % len(in_multiple_fovs))

        # find the boundaries of the clipped cells
        _list = self.collate_borders_par(in_multiple_fovs)
        df_2 = pd.DataFrame(_list).astype({"label": int})

        # Both clipped and unclipped in a dataframe
        res = pd.concat([df_1, df_2])

        set_diff = set(cell_props.label) - set(res.label.values)
        if set_diff:
            unresolved_labels = pd.DataFrame({'label': list(set_diff), 'coords': np.nan * np.ones(len(set_diff))});
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
        return out

    def find_offset(self, fov_ids):
        sanity_check = np.array([self.fovs[d]['fov_id'] == d for d in fov_ids])
        assert np.all(sanity_check)
        offset_x = min([self.fovs[d]['fov_offset_x'] for d in fov_ids])
        offset_y = min([self.fovs[d]['fov_offset_y'] for d in fov_ids])
        return offset_x, offset_y

    def obj_outline(self, fov, cell_props):
        logger.info('Getting cell boundaries for cells in fov: %d' % fov['fov_id'])
        label_image = fov['label_image'].toarray()
        offset_x = fov['fov_offset_x']
        offset_y = fov['fov_offset_y']
        clipped_cells = cell_props[cell_props.is_clipped].label.values

        df = extract_borders_par(label_image, offset_x, offset_y, clipped_cells)
        return df

    def poly_vertices(self, my_list):
        '''
        Removes points that exist on the same line. Ie it keeps on the vertices along the perimeter of the polygon
        :param my_list: list of tuples. Has to be sorted, either clockwise or counterclock wise. I do not check if this
        condition is met.
        :return:
        '''
        for i in reversed(range(1, len(my_list) - 1)):
            a = my_list[i - 1]
            x = my_list[i]
            b = my_list[i + 1]
            if (x[1] - a[1]) * (b[0] - a[0]) == (b[1] - a[1]) * (x[0] - a[0]):
                del my_list[i]
        return my_list

    def sort_ccw(self, l):
        '''
        Sort the list of coords counter-clockwise.
        BE CAREFUL!!! The origin is assumed to be the top-left corner. Remove the 'reverse=True' flag if your origin is
        at the bottom left corner
        TO CHECK: WHAT HAPPENS IF YOU SORT AN ALREADY SORTED ARRAY?
        In general I need to test this code a bit more
        :param l:
        :return:
        '''
        # remove duplicates from the list
        l = list(set(l))
        yAvg = sum(t[1] for t in l) / len(l)
        xAvg = sum(t[0] for t in l) / len(l)

        def algo(t):
            return (math.atan2(t[1] - yAvg, t[0] - xAvg) + 2 * math.pi) % (2 * math.pi)

        # sort the list counterclockwise
        l.sort(key=algo, reverse=True)  # <--- REMOVE THE 'reverse=True' flag if your origin is the bottom-left corner

        # last point same as the first one
        l.extend([l[0]])
        # logger.info('Coords sorted counter clockwise')
        return l


if __name__ == "__main__":
    fovs_across = 11
    fovs_down = 14

    fovs = Fov(fovs_across, fovs_down, config.PREPROCESSOR)
    stage = Stage(fovs)
    # stage_spots = Stage_spots(fovs)

    stage.merge_cells()
    stage.global_labels()
    stage.cell_props = stage.calc_props()

    logger.info('Spot labelling started')
    res_list = []
    for i, fov in enumerate(stage.fovs):
        assert i == fov['fov_id'], 'The list is miss-aligned'
        # Assign the parent cell and its coordinates
        res = stage.assign_spot_parent(fov)
        fov['spots'] = res
        res['fov_id'] = fov['fov_id']
        res_list.append(res)
    logger.info('Spot labelling done')

    stage.cell_props['cell_id'] = stage.assign_cell_id()

    stage.spots = pd.concat(res_list).astype({"label": int})

    # Save now the data on the filesystem
    # stage.writer()

    print('Done!')