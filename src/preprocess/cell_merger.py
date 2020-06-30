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


# def _get_connected_labels(labels):
#     connected_labels = []
#     label_set = set(labels)
#     for label in list(label_set):
#         temp = [i for i, j in enumerate(labels) if j == label]
#         if len(temp) > 1:
#             connected_labels.append(np.array(temp))
#     return np.array(connected_labels)
#
#
# def _get_connected_labels_fast(labels):
#     positions = {}
#     for index, label in enumerate(labels):
#         positions.setdefault(label, []).append(index)
#     return np.array([np.array(indices) for indices in positions.values() if len(indices) > 1])


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



def load_spots(fov_id, scaling_factor, my_dir):
    '''
    :param fov: array of the fov_ids that the block consists of
    :return:
    '''
    # fov_id = fov['fov_id']
    # my_dir = self.my_config['STARFISH_SPOTS']  # That has be kept in a more tidy manner. It is a bit of a mess!

    # fov_csv = [f"fov_{int(fov_id):03d}.csv"]
    fov_csv = ['spots_%d.csv' % fov_id]
    filename = [os.path.join(my_dir, f) for f in os.listdir(my_dir) if f in fov_csv]

    if filename:
        f = filename[0]
    else:
        f = None

    logger.info('reading: %s' % f)
    spots = pd.read_csv(f).dropna()

    spots = spots[['Gene', 'x', 'y']]
    spots = spots[~spots.duplicated()]
    gene_name, idx = np.unique(spots.Gene.values, return_inverse=True)
    spots['gene_id'] = idx  # this is effectively the gene id

    spots['x'] = spots.x * scaling_factor
    spots['y'] = spots.y * scaling_factor
    spots = spots.sort_values(['x', 'y'], ascending=[True, True]) \
        .reset_index(drop=True)  # <-- DO NOT FORGET TO RESET THE INDEX

    return spots


def update_spots(fov, spots):
    fov['spots'] = spots
    return fov


class Stage(object):
    def __init__(self, fovs_obj, cellmaps):
        self._fovs_obj = fovs_obj
        self.my_counter = itertools.count()
        self.merge_register = defaultdict(list) # keeps the new created labels (as keys) and the fovs where these labels appear (as values)
        self.cell_props = None
        self.spots = None
        self.cellmaps = cellmaps

        # file_list = [self.make_path(d) for d in range(len(self.fovs))]
        # start = time.time()
        # for i, d in enumerate(self._fovs_obj.fovs[:20]):
        #     d['label_image_orig'] = self.load_label_image(i)
        #     d['spots_orig'] = load_spots(i, self.scaling_factor, self.my_config['MATLAB_SPOTS'])
        # logger.info('Finished Single thread in %s secs' % (time.time()-start))

        start = time.time()
        res = self.my_multithread(list(range(len(self._fovs_obj.fovs))))
        for i, d in enumerate(self._fovs_obj.fovs):
            d['label_image'] = res[i][0]
            d['spots'] = res[i][1]
        logger.info('Finished Multi thread in %s secs' % (time.time() - start))
        print('after multithread')

    @property
    def fovs(self):
        return self._fovs_obj.fovs

    @property
    def fovs_across(self):
        return self._fovs_obj.fovs_across

    @property
    def fovs_down(self):
        return self._fovs_obj.fovs_down

    @property
    def fov_shape(self):
        return self._fovs_obj.fov_shape

    @property
    def scaling_factor(self):
        return self._fovs_obj.scaling_factor

    @property
    def my_config(self):
        return self._fovs_obj.my_config

    def load(self, i):
        if self.cellmaps is not None:
            label_image = coo_matrix(self.cellmaps[i])
        else:
            label_image = self.load_label_image(i)
        spots = load_spots(i, self.scaling_factor, self.my_config['MATLAB_SPOTS'])
        return label_image, spots

    def my_multithread(self, ids):
        n = max(1, cpu_count() - 1)
        pool = ProcessPool(n)
        results = pool.map(self.load, ids)
        pool.close()
        pool.join()
        return results

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
        # img_obj = fov['image_objects']
        # img_obj = unlocalise_coords(img_obj, block)
        # img_obj['x_cell'] = img_obj.x
        # img_obj['y_cell'] = img_obj.y
        # img_obj = img_obj.drop(columns=['x', 'y'])
        #
        # # img_obj.index = img_obj.index.astype('Int64')
        # # spots['label'] = spots.label.astype('Int64')
        # spots = spots.merge(img_obj, left_on='label', right_index=True, how='left')

        out['label'] = out.label.astype('Int64')
        return out

    def localise_coords(self, spots, fov):
        '''
        convert spots to local
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

    def get_tile_coords_fake(self, f):
        out = []
        out.append(([0, 10], [0, 10], 1.0, np.array([11, 11])))
        out.append(([10, 20], [0, 10], 1.0, np.array([11, 11])))
        out.append(([0, 10], [10, 20], 1.0, np.array([11, 11])))
        out.append(([10, 20], [10, 20], 1.0, np.array([11, 11])))
        return out[f]

    # def global_coords(self, coords):
    #     x_range = coords[0]
    #     y_range = coords[1]
    #     k = coords[2]
    #     x_global = [d * k for d in x_range]
    #     y_global = [d * k for d in y_range]
    #     return x_global, y_global, k

    # def get_dir(cfg, fov_id):
    #     return os.path.join(cfg, 'fov_' + str(fov_id))

    def get_tile_coords(self, fov):
        '''
        calc the tile coordinates and the scaling factor of the tile coordinates. Multiply by the scaling factor
        to get the actual global pixel coordinates. For example if:

            xc = [650.0, 975.0]
            yc = [0.0, 325.0]
            scaling_factor = 6.153846153846154

            then the pixel coordinates are
            xc = [4000.0, 6000.0]
            yc = [0.0, 2000.0]

        :param fov:
        :return:
        '''

        nuclei_dir = r"https://raw.githubusercontent.com/acycliq/starfish_ISS_h_brain_03/master/main_files/"
        jsn = f"nuclei-fov_{int(fov):03d}.json"

        request = urllib.request.Request(urllib.parse.urljoin(nuclei_dir, jsn))
        base64string = base64.b64encode(bytes('%s:%s' % (credentials.USER, credentials.TOKEN), 'ascii'))
        request.add_header("Authorization", "Basic %s" % base64string.decode('utf-8'))
        with urllib.request.urlopen(request) as result:
            logger.info('loading..')
            data = json.loads(result.read().decode())
            logger.info('loaded..')
            # data = json.load(result) # i think that also works, dont know whats best!
            xc = data['tiles'][0]['coordinates']['xc']
            yc = data['tiles'][0]['coordinates']['yc']
            default_tile_shape = np.array(data['default_tile_shape'])
            tile_shape = np.array(data['tiles'][0]['tile_shape'])
            scale_x = data['default_tile_shape'][0] / np.diff(data['tiles'][0]['coordinates']['xc'])[0]
            scale_y = data['default_tile_shape'][1] / np.diff(data['tiles'][0]['coordinates']['yc'])[0]
            assert (scale_y == scale_x)
            assert (np.all(default_tile_shape == tile_shape))
            scaling_factor = scale_x

        return xc, yc, scaling_factor, tile_shape

    def load_label_image(self, fov_id):
        '''
        Reads a label_image and returns a dataframe with the objects found on it.
        The array has index the label and columns: fov_id, area and the cell centroids x, y
        :param fov:
        :return:
        '''
        # fov = self.fovs[fov_id]
        # fov_id = fov['fov_id']
        label_img = self.read_from_disk(fov_id)
        return coo_matrix(label_img)

    def read_from_disk(self, fov_id):
        full_path = self.make_path(fov_id)
        # print('reading %s' % full_path)
        out = np.genfromtxt(full_path, delimiter=',')
        logger.info('reading %s', full_path)
        return out

    def make_path(self, fov_id):
        # header_flag is None is no header exist
        # x could be 'label_image'
        fov_dir = get_dir(self.my_config, fov_id)
        full_path = os.path.join(fov_dir, 'label_image', 'label_image_fov_' + str(fov_id) + '.csv')
        return full_path

        # fov_id_str = 'fov_%d' % n
        # fName_str = 'label_image_fov_%d.csv' % n
        # path_str = os.path.join('.', 'demo_data', 'human_full', 'fovs', fov_id_str, 'raw', 'label_image', fName_str)
        # return path_str

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
        for fov in self.fovs:
            logger.info('\n')
            logger.info('Doing fov %i' % fov['fov_id'])
            self.merge(fov)
        logger.info('Relabelling finished')
        logger.info('\n')

    def merge(self, fov):
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

        # maxab = max(temp_a.data.max(), temp_b.data.max())
        # my_list = list(range(maxab, maxab + len(border)))
        for d in matched_labels:
            new_label = self._new_label(d)
            for x in d['a']:
                temp_a.data[temp_a.data == x] = new_label
                self.update_register(adjc_fov['fov_id'], new_label, x)
                # logger.info('fov_%d: label %d ---> label %d' % (adjc_fov['fov_id'], x, new_label))
                # sleep(0.1)

            for x in d['b']:
                temp_b.data[temp_b.data == x] = new_label
                self.update_register(fov['fov_id'], new_label, x)
                # logger.info('fov_%d: label %d ---> label %d' % (fov['fov_id'], x, new_label))
                # sleep(0.1)

        # print('fov_id: (%d, %d), %s' % (adjc_fov['fov_id'], fov['fov_id'], dlist))
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
        '''
        Reads a label_image and returns a dataframe with the objects found on it.
        The array has index the label and columns: fov_id, area and the cell centroids x, y
        :param fov:
        :return:
        '''
        fov_id = fov['fov_id']
        label_image = fov['label_image'].toarray().astype(np.int32)

        props = skmeas.regionprops(label_image)
        props_arr = [[int(p.label), int(fov_id), int(p.area), p.centroid[1], p.centroid[0]] for i, p in
                     enumerate(props)]
        props_df = pd.DataFrame(props_arr, columns=['label', 'fov_id', 'area', 'x_local', 'y_local'])

        # stick the df to the dict
        fov['image_objects'] = props_df  # <--- Not sure If I should do this or just return the df and add the dictionary filed at the calliing function
        return fov

    def global_labels_par(self):
        results = self._global_labels_par()

        # finally flip the sign on the dict
        self.flip_sign()
        return results

    def _global_labels_par(self):
        n = max(1, cpu_count() - 1)
        pool = ThreadPool(n)
        results = pool.map(self.global_labels_helper, self.fovs)
        pool.close()
        pool.join()

        return results

    def global_labels(self):
        '''
        Same as "global_labels_par()" but without parallelism
        :return:
        '''
        for fov in self.fovs:
            # logger.info('fov:%d, Setting global labels' % fov['fov_id'])
            data = fov['label_image'].data
            labels = np.unique(data[data > 0])
            label_map = {d: self.label_generator() for d in labels}
            fov['label_image'].data = np.array([label_map[d] if d > 0 else d for d in data])
            logger.info('fov: %d: label map is: %s' % (fov['fov_id'], label_map))

            clipped_cell_labels = {k: v for k, v in self.merge_register.items() if fov['fov_id'] in v}.keys()
            assert np.all([d in fov['label_image'].data for d in clipped_cell_labels]), \
                "A label that the cell register says should exist in the fov %d doesnt seem to be verified by the label_image" % \
                fov['fov_id']

            fov['label_image'].data = -1 * fov['label_image'].data
            if np.any(fov['label_image'].data):
                logger.info('fov:%d, Global labels are set.' % fov['fov_id'])
            else:
                logger.info('fov:%d empty, nothing to do here' % fov['fov_id'])

            self.image_objects(fov)

        # finally flip the sign on the dict
        self.flip_sign()


    def global_labels_helper(self, fov):
        logger.info('fod_id is: %d' % fov['fov_id'])
        data = fov['label_image'].data
        labels = np.unique(data[data > 0])
        label_map = {d: self.label_generator() for d in labels}
        fov['label_image'].data = np.array([label_map[d] if d > 0 else d for d in data])

        # Sanity check. merge register keeps the cells that are clipped, the cells that span over more than a single fov
        # Get the labels from the merge register for this particular fov
        clipped_cell_labels = {k: v for k, v in self.merge_register.items() if fov['fov_id'] in v}.keys()
        assert np.all([d in fov['label_image'].data for d in clipped_cell_labels]), \
            "A label that the cell register says should exist in the fov %d doesnt seem to be verified by the label_image" % fov['fov_id']

        # Flip now the sign of the labels
        fov['label_image'].data = -1 * fov['label_image'].data
        if np.any(fov['label_image'].data):
            logger.info('fov:%d, Global labels are set.' % fov['fov_id'])
        else:
            logger.info('fov:%d empty, nothing to do here' % fov['fov_id'])

        self.image_objects(fov)

        return True

    def flip_sign(self):
        # need to flip the sign of the keys in the merge_register dict
        # Maybe I should add some sanity checks here.
        # For example I can do the flipping inside the loop of the calling function and while you loop get the keys (ie the labels)
        # for the particular fov and then flip the label sign (if it negative because it may have been flipped already if the label
        # exists in more that one fov

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
        '''
        Reads a label_image and returns a dataframe with the objects found on it.
        The array has index the label and columns: fov_id, area and the cell centroids x, y
        :param fov:
        :return:
        '''

        logger.info('merged_props_dict2 starts')
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

            logger.info('calculating centroid and area for the merged cell with label %s' % labels)
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

                logger.info('label %d spans across fov_ids %s' % (p.label, fov_ids))
                # if (p.label == 1450):
                #     stop_here = 1
                # logger.info('label %d centroid is (x,y): %s, area: %4.2f' % (p.label, p.centroid[::-1], p.area))
                # logger.info('label %d centroid local coords are (x, y): %s ' % (p.label, centroid_coords[::-1]))
                # logger.info('label %d centroid falls within the fov with id: %d' % (p.label, centroid_fov_id))

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
                logger.info('merged_props_dict ends')

        # # first, do the merge cells
        # # check the register to see which cells are merged
        # merged_props_dict = {'label': [], 'fov_id': [], 'area': [], 'x_local': [], 'y_local': []}
        # for t in self.merge_register.items():
        #     label = t[0]
        #     fov_ids = sorted(t[1])
        #
        #     logger.info('calculating centroid and area for the merged cell with label %d' % label)
        #     if len(fov_ids) > 1:
        #         label_image = self.collate_arrays(fov_ids)
        #     else:
        #         label_image = None  # you shouldnt get here. Since we merge cells we need at least two fovs
        #
        #     props = skmeas.regionprops(label_image.astype(np.int32))
        #     p = [d for d in props if d.label == label][0]
        #     centroid_fov_id, centroid_coords = self.locate_fov(p.centroid, fov_ids)
        #
        #     logger.info('label %d spans across fov_ids %s' % (p.label, fov_ids))
        #     # logger.info('label %d centroid is (x,y): %s, area: %4.2f' % (p.label, p.centroid[::-1], p.area))
        #     # logger.info('label %d centroid local coords are (x, y): %s ' % (p.label, centroid_coords[::-1]))
        #     # logger.info('label %d centroid falls within the fov with id: %d' % (p.label, centroid_fov_id))
        #
        #     # sanity check
        #     filtered_fovs = list(filter(lambda d: d['fov_id'] in fov_ids, self.fovs))
        #     mask = [d['image_objects'].label == label for d in filtered_fovs]
        #     total_area = sum([d['image_objects'].area[mask[i]].values[0] for i, d in enumerate(filtered_fovs)])
        #     # logger.info('total sum for the sanity check is: %d' % total_area)
        #     assert int(total_area) == int(p.area)
        #
        #     merged_props_dict['label'].append(p.label)
        #     merged_props_dict['fov_id'].append(centroid_fov_id)
        #     merged_props_dict['area'].append(p.area)
        #     merged_props_dict['x_local'].append(centroid_coords[1])
        #     merged_props_dict['y_local'].append(centroid_coords[0])

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

        # logger.info('idx1 start')
        # idx1 = np.where(_x[:, None] == _y[None, :])[1]
        # logger.info('idx1 end')
        #
        # logger.info('idx2 start')
        # idx2 = [i for d in _x for i, el in enumerate(_y) if el == d]
        # logger.info('idx2 end')

        # logger.info('idx3 start')
        idx = self.my_remap(_x, _y)
        # logger.info('idx3 end')

        # assert np.all(idx1 == idx2)
        # assert np.all(idx1 == idx3)
        return label_map[idx, 1]
        # print('hello')

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

        # coords = {}
        # coords['x'] = coord_row
        # coords['y'] = coord_col
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
        # a = np.arange(20).reshape(4, 5)
        mask = np.isin(a, d)
        return np.where(mask, a, np.nan)[mask.any(axis=1)][:, mask.any(axis=0)]

    def assign_cell_id(self):
        # finally save the cell props to a csv file
        # Need to add the cell_id. This should be made redundant. The label can be used instead.
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
        cell_props[cells_headers].to_csv('expanded_cells_david.csv', index=False)

        # 2. save the cell coords
        coords_headers = ['cell_id', 'label', 'coords']
        cell_props[coords_headers].to_json('cell_coords_david.json', orient='records')

        # 3. save the spots
        spots_df = self.spots.copy()
        spots_df['target'] = spots_df.Gene
        spots_df['x_global'] = spots_df.x
        spots_df['y_global'] = spots_df.y
        spots_df['fov_id'] = spots_df.fov_id
        spots_df['x_cell'] = spots_df.x_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)
        spots_df['y_cell'] = spots_df.y_cell.fillna(-1).astype(int).astype('str').replace('-1', np.nan)

        spots_headers = ['x_global', 'y_global', 'fov_id', 'label', 'target', 'x_cell', 'y_cell']
        spots_df[spots_headers].to_csv('spots_david.csv', index=False)
        logger.info('Total number of collected spots: %d' % spots_df.shape[0])

        # # 3b. Save now the spots seperately for each fov
        # for fov in self.fovs:
        #     fov_id = fov['fov_id']
        #     fov_dir = get_dir(self.my_config, fov_id)
        #     full_path = os.path.join(fov_dir, 'spots', 'spots_fov_' + str(fov_id) + '.csv')
        #
        #     if not os.path.exists(os.path.dirname(full_path)):
        #         os.makedirs(os.path.dirname(full_path))
        #     df = fov['spots']
        #     df['target'] = df.Gene
        #     df['x_global'] = df.x
        #     df['y_global'] = df.y
        #     df['fov_id'] = df.fov_id
        #     # target, x_global, y_global, Expt, label, x_cell, y_cell
        #     df = df[['x_global', 'y_global', 'fov_id', 'gene_id', 'label', 'target', 'x_cell', 'y_cell']]
        #     df.to_csv(full_path, index=False)

    def cell_boundaries(self, cell_props):
        '''
        calculate the outlines of the cells
        :return:
        '''

        # loop over the self.cell_props
        res_list = []
        for fov in self.fovs:
            if np.any(fov['label_image'].data):
                df = self._obj_outline(fov, cell_props)
                res_list.append(df)
            else:
                logger.info('fov:%d empty, No cells to draw boundaries were found' % fov['fov_id'])
                # # label_image = fov['label_image']
        _df = pd.concat(res_list).astype({"label": int})

        # check if the cell spans across more than one fov
        # dict_keys = np.array(list(self.merge_register.keys()))
        # idx = self.my_remap(_df.label.values, dict_keys).mask
        # df_1 = _df[idx]  # <--- Dataframe the keeps boundaries of the cells which are not clipped by the fov
        df_1 = _df.iloc[np.isin(_df.label, cell_props[~cell_props.is_clipped].label)]  # <--- Dataframe the keeps boundaries of the cells which are not clipped by the fov

        in_multiple_fovs = sorted(cell_props[cell_props.is_clipped].label.values)
        logger.info('There are %d cells whose boundaries span across multiple fovs' % len(in_multiple_fovs))
        # _list = []
        # logger.info('Single Thread starts')
        # for label in in_multiple_fovs:
        #     logger.info('label: %d. Finding the cell boundaries' % label)
        #     _d = self.merge_register
        #     label_image = self.collate_arrays(_d[label])
        #     offset_x, offset_y = self.find_offset(_d[label])
        #     out = self._obj_outline_helper(label_image, offset_x, offset_y)
        #     # keep only the outline relevant to the label that spans across several fovs
        #     out = out[out.label == label]
        #     _list.append(out)
        # df_2 = pd.concat(_list).astype({"label": int})
        # logger.info('Single Thread ends')
        # logger.info('Multi Thread starts')
        _list = self.collate_borders_par(in_multiple_fovs)
        df_2 = pd.concat(_list).astype({"label": int})
        # logger.info('Multi Thread ends')
        res = pd.concat([df_1, df_2])

        # sort now the coords countrerclockwise
        _coords = res.coords.map(self.sort_ccw)

        # remove now redundant coordinates
        _coords = _coords.map(self.poly_vertices)

        # replace now the coords column with the clean one
        res['coords'] = _coords
        logger.info('Coords sorted counter clockwise')

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
        logger.info('label: %d. Finding the cell boundaries' % label)
        label_image = self.collate_arrays(self.merge_register[label])
        offset_x, offset_y = self.find_offset(self.merge_register[label])
        out = self._obj_outline_helper(label_image, offset_x, offset_y)
        # keep only the outline relevant to the label that spans across several fovs
        out = out[out.label == label]
        return out

    def find_offset(self, fov_ids):
        sanity_check = np.array([self.fovs[d]['fov_id'] == d for d in fov_ids])
        assert np.all(sanity_check)
        offset_x = min([self.fovs[d]['fov_offset_x'] for d in fov_ids])
        offset_y = min([self.fovs[d]['fov_offset_y'] for d in fov_ids])
        return offset_x, offset_y

    def _obj_outline(self, fov, cell_props):
        label_image = fov['label_image'].toarray()
        offset_x = fov['fov_offset_x']
        offset_y = fov['fov_offset_y']
        df = self._obj_outline_helper(label_image.copy(), offset_x, offset_y)

        # Filter now the df and keep only the labels that also exist in cell_props for that fov and are
        # not clipped.
        # Also, If a label exists in cell_props (filtered as now explained) and not in df then raise
        # a warning that the cell has not outline
        mask = (cell_props.fov_id == fov['fov_id']) & (~cell_props.is_clipped)
        cell_props_masked = cell_props[mask]
        df_masked = df.iloc[np.isin(df.label, cell_props_masked.label)]

        set_diff = set(cell_props_masked.label.values) - set(df_masked.label)
        loop_depth = 3
        df_list = []
        df_list.append(df)

        recurse_depth = 1
        while (set_diff and recurse_depth < 5):
            logger.info('Doing another pass because these nested cells were found %s', set_diff)
            logger.info('Pass Num: %d' % recurse_depth)
            df_2 = self._obj_outline_helper(label_image.copy(), offset_x, offset_y, set_diff)
            df_list.append(df_2)
            df = pd.concat(df_list).astype({"label": int})
            df_masked = df.iloc[np.isin(df.label, cell_props_masked.label)]
            set_diff = set(cell_props_masked.label.values) - set(df_masked.label)
            recurse_depth = recurse_depth + 1


        if set_diff:
            logger.info('Couldnt derive the boundaries for cells with label: %s' % list(set_diff))

        return df


        # if np.any(fov['label_image'].data):
        #     set_diff = set(fov['label_image'].data) - set(df.label.values)
        #     if set_diff:
        #         logger.info('Couldnt derive the boundaries for cells with label: %s' % list(set_diff))




    def _obj_outline_helper(self, label_image, offset_x, offset_y, keep_list=None):
        if keep_list:
            for label in set(label_image.flatten()) - set(keep_list):
                label_image[label_image == label] = 0
        mask = ndimage.binary_erosion(label_image)
        label_image[mask] = 0
        c = coo_matrix(label_image)
        # transform to global coords
        c_col = c.col.astype(int) + int(offset_x)
        c_row = c.row.astype(int) + int(offset_y)
        if c.data.size > 0:
            my_df = pd.DataFrame({'coords': list(zip(c_col, c_row)), 'label': c.data})
            my_df = my_df.groupby('label')['coords'].apply(lambda group_series: group_series.tolist()).reset_index()
            my_df = my_df.astype({"label": int})
        else:
            logger.info('array empty')
            my_df = pd.DataFrame()
        return my_df


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