import numpy as np
import pandas as pd
import json
import numexpr as ne
import numba as nb
import scipy
import xarray as xr
import scipy.io as spio
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix, hstack
import os
import config
import urllib, base64
from itertools import accumulate
# import credentials
import gc
import logging

dir_path = os.path.dirname(os.path.realpath(__file__))
CONFIG_FILE = dir_path + '/config.yml'


logger = logging.getLogger()


def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects

    from: `StackOverflow <http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries>`_
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


def label_spot(a, idx):
    '''
    Given an array a (image_array) and
    :param a: An array of size numPixelsY-by-numPixelsX specifying that element (i,j) belongs to
                cell a[i,j]. Note that array a is 1-based, ie if pixel (i,j) is outside a cell then
                a[i,j] = 0.
    :param idx: An array of size 2-by-N of the pixels coordinates of spot idx[k], k=1...N
    :return:
    a = np.array([  [4,0,1],
                    [2,0,0],
                    [0,1,0]])

    idx = np.array([[0,0],
                    [2, 1],
                    [1,2],
                    [1,3]])

    IndexArrayNan(a, idx.T) = [4., 1., 0., nan]
    which means that:
            spot with coords [0,0] belongs to cell 4
            spot with coords [2,0] belongs to cell 1
            spot with coords [1,2] belongs to 0 (ie no assigned cell)
            spot with coords [1,3] is outside the bounds and assigned to nan

    '''
    assert isinstance(idx[0], np.ndarray), "Array 'idx' must be an array of arrays."
    idx = idx.astype(np.int64)
    out = np.array([])
    dim = np.ones(idx.shape[0], dtype=int)
    dim[:len(a.shape)] = a.shape

    # output array
    out = np.nan * np.ones(idx.shape[1], dtype=int)

    # find the ones within bounds:
    is_within = np.all(idx.T <= dim-1, axis=1)

    # also keep only non-negative ones
    is_positive = np.all(idx.T >= 0, axis=1)

    # filter array`
    arr = idx[:, is_within & is_positive]
    flat_idx = np.ravel_multi_index(arr, dims=dim, order='C')
    out[is_within & is_positive] = a.ravel()[flat_idx]

    print('in label_spot')

    return out


def gammaExpectation(rho, beta):
    '''
    :param r:
    :param b:
    :return: Expectetation of a rv X following a Gamma(r,b) distribution with pdf
    f(x;\alpha ,\beta )= \frac{\beta^r}{\Gamma(r)} x^{r-1}e^{-\beta x}
    '''

    # sanity check
    assert (np.all(rho.coords['cell_id'].data == beta.coords['cell_id'])), 'rho and beta are not aligned'
    assert (np.all(rho.coords['gene_name'].data == beta.coords['gene_name'])), 'rho and beta are not aligned'
    r = rho.data[:, :, None]
    b = beta.data
    gamma = np.empty(b.shape)
    ne.evaluate('r/b', out=gamma)
    out = xr.DataArray(gamma.astype(np.float64),
                       coords={'cell_id': beta.cell_id.values,
                               'gene_name': beta.gene_name.values,
                               'class_name': beta.class_name.values},
                       dims=['cell_id', 'gene_name', 'class_name']
                       )
    del gamma
    del r
    del b
    gc.collect()
    del gc.garbage[:]
    return out


def logGammaExpectation(rho, beta):
    '''
    :param r:
    :param b:
    :return: Expectetation of a rv log(X) where X follows a Gamma(r,b) distribution
    '''
    # start = time.time()
    # out = scipy.special.psi(r) - np.log(b)
    # end = time.time()
    # print('time in logGammaExpectation:', end - start)
    # start = time.time()

    # sanity check
    assert (np.all(rho.coords['cell_id'].data == beta.coords['cell_id'])), 'rho and beta are not aligned'
    assert (np.all(rho.coords['gene_name'].data == beta.coords['gene_name'])), 'rho and beta are not aligned'
    r = rho.data[:, :, None]
    b = beta.data
    logb = np.empty(b.shape)
    ne.evaluate("log(b)", out=logb)
    log_gamma = scipy.special.psi(r) - logb

    out = xr.DataArray(log_gamma.astype(np.float64),
                       coords={'cell_id': beta.cell_id.values,
                               'gene_name': beta.gene_name.values,
                               'class_name': beta.class_name.values},
                       dims=['cell_id', 'gene_name', 'class_name']
                       )
    # end = time.time()
    # print('time in logGammaExpectation ne:', end - start)

    del logb
    del log_gamma
    del r
    del b
    gc.collect()
    del gc.garbage[:]

    return out


def negBinLoglik(da_x, r, da_p):
    '''
    Negative Binomial loglikehood
    :param x:
    :param r:
    :param p:
    :return:
    '''

    # sanity check
    assert (np.all(da_x.coords['cell_id'].data == da_p.coords['cell_id'])), 'gene counts and beta probabilities are not aligned'
    assert (np.all(da_x.coords['gene_name'].data == da_p.coords['gene_name'])), 'gene counts and beta probabilities are not aligned'

    contr = np.zeros(da_p.shape)
    x = da_x.data[:, :, None]
    p = da_p.data
    # start = time.time()
    # out = x * np.log(p, where=x.astype(bool)) + r * np.log(1-p)
    # end = time.time()
    # print('time in negBinLoglik:', end - start)
    # start = time.time()
    ne.evaluate("x * log(p) + r * log(1 - p)", out=contr)
    # end = time.time()
    # print('time in negBinLoglik - ne:', end - start)

    out = xr.DataArray(contr,
                       coords={'cell_id': da_p.cell_id.values,
                               'gene_name': da_p.gene_name.values,
                               'class_name': da_p.class_name.values},
                       dims=['cell_id', 'gene_name', 'class_name']
                       )

    return out


def _negBinLoglik(x, r, p):
    '''
    Negative Binomial loglikehood
    :param x:
    :param r:
    :param p:
    :return:
    '''
    out=np.zeros(p.shape)
    # start = time.time()
    # out = x * np.log(p, where=x.astype(bool)) + r * np.log(1-p)
    # end = time.time()
    # print('time in negBinLoglik:', end - start)
    # start = time.time()
    ne.evaluate("x * log(p) + r * log(1 - p)", out=out)
    # end = time.time()
    # print('time in negBinLoglik - ne:', end - start)
    return out

@nb.njit(parallel=True, fastmath=True)
def nb_negBinLoglik(x, r, p):
    '''
    Negative Binomial loglikehood
    :param x:
    :param r:
    :param p:
    :return:
    '''
    out = np.empty(p.shape,p.dtype)

    for i in nb.prange(p.shape[0]):
        for j in range(p.shape[1]):
            if x[i, j, 0] != 0.:
                x_ = x[i, j, 0]
                for k in range(p.shape[2]):
                    out[i, j, k] = x_ * np.log(p[i, j, k]) + r * np.log(1.-p[i, j, k])
            else:
                for k in range(p.shape[2]):
                    out[i, j, k] = r * np.log(1.-p[i, j, k])

    return out


def softmax(x):
    '''
    https://stackoverflow.com/questions/34968722/how-to-implement-the-softmax-function-in-python
    :param x:
    :return:
    '''
    """Compute softmax values for each sets of scores in x."""
    assert 'class_name' in x.dims, 'There is not dimension "class_name" '
    temp = xr.ufuncs.exp(x)
    return temp / temp.sum('class_name')


def softmax2(x):
    '''
    https://stackoverflow.com/questions/34968722/how-to-implement-the-softmax-function-in-python
    :param x:
    :return:
    '''
    """Compute softmax values for each sets of scores in x."""
    return np.exp(x) / np.sum(np.exp(x), axis=1)[:, None]


def softmax_nolan(X, theta = 1.0, axis = None):
    """
    From https://nolanbconaway.github.io/blog/2017/softmax-numpy
    Compute the softmax of each element along an axis of X.

    Parameters
    ----------
    X: ND-Array. Probably should be floats.
    theta (optional): float parameter, used as a multiplier
        prior to exponentiation. Default = 1.0
    axis (optional): axis to compute values along. Default is the
        first non-singleton axis.

    Returns an array the same size as X. The result will sum to 1
    along the specified axis.
    """

    # make X at least 2d
    y = np.atleast_2d(X)

    # find axis
    if axis is None:
        axis = next(j[0] for j in enumerate(y.shape) if j[1] > 1)

    # multiply y against the theta parameter,
    y = y * float(theta)

    # subtract the max for numerical stability
    y = y - np.expand_dims(np.max(y, axis = axis), axis)

    # exponentiate y
    y = np.exp(y)

    # take the sum along the specified axis
    ax_sum = np.expand_dims(np.sum(y, axis = axis), axis)

    # finally: divide elementwise
    p = y / ax_sum

    # flatten if X was 1D
    if len(X.shape) == 1: p = p.flatten()

    return p


def isConverged(spots, p0, tol):
    p1 = spots.call.cell_prob.values
    if p0 is None:
        p0 = np.zeros(p1.shape)
    delta = np.max(np.abs(p1 - p0))
    converged = (delta < tol)
    return converged, delta


def bi2(X, dims, *args):
    nK = dims[1]
    inds = []
    if len(args) == 2:
        # print('im in!')
        args = (*args, np.arange(0, nK))

    temp = np.zeros(dims).astype(int)
    for i in range(len(args)):
        inds.append(temp + args[i])

    inds = tuple(inds)
    # out = X[inds]
    return X[inds]


def rename_header(rootdir):
    '''
    Convenience function to rename the headef of a csv file.
    I used this for the expanded_cell csv flatfiles.
    :param rootdir:
    :return:
    '''
    header = 'cell_id,metadata_position,area,x,y\n'
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            f=open(os.path.join(rootdir, file),'r') #use the absolute URL of the file
            lines = f.readlines()
            f.close()
            f = open(file,'w') #universal mode can only be used with 'r' mode
            for line in lines:
                if line == header:
                    line = 'cell_id,fov_id,area,row,col\n'
                f.write(line)
            f.close()



def collate_expanded_cells():
    '''
    reads cell segmentation files (csv) for each fov and then concatenates and saves in a csv
    the cells' centroid coordinates, their area and which field of view the corresponding cell
    appears in
    '''

    logger = logging.getLogger()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s:%(levelname)s:%(message)s"
    )

    ROOT_DIR = config.ROOT_DIR
    my_dir = os.path.join(ROOT_DIR, 'demo_data', 'human', 'expanded_cells')
    filenames = [os.path.join(my_dir, f) for f in os.listdir(my_dir) if f.endswith('.' + 'csv')]

    logger.info('joining csvs...')
    df = pd.concat([pd.read_csv(f) for f in filenames], sort=False)
    fov_id = ['fov_%03d' % x for x in df.fov_id]
    df['cell_id'] = np.arange(df.shape[0])  # overwrite cell_id  because it gets reset every time a new csv is parsed.
    df['fov_id'] = fov_id
    # print(df.head())
    # print(df.shape)

    nuclei_dir = r"https://raw.githubusercontent.com/acycliq/starfish_ISS_h_brain_03/master/main_files/"
    my_list = []
    target_url = urllib.parse.urljoin(nuclei_dir, 'nuclei.json')

    with urllib.request.urlopen(target_url) as url:
        d = json.loads(url.read().decode())
        my_data = d['contents']
        for fov, jsn in my_data.items():
            with urllib.request.urlopen(urllib.parse.urljoin(nuclei_dir, jsn)) as f:
                data = json.loads(f.read().decode())
                # data = json.load(f)
                x_0 = data['tiles'][0]['coordinates']['xc'][0]
                y_0 = data['tiles'][0]['coordinates']['yc'][0]
                scale_x = data['default_tile_shape'][0] / np.diff(data['tiles'][0]['coordinates']['xc'])[0]
                scale_y = data['default_tile_shape'][1] / np.diff(data['tiles'][0]['coordinates']['yc'])[0]
                my_list.append([fov, x_0, y_0, scale_x, scale_y])

    tile_origin = pd.DataFrame(my_list, columns=['fov_id', 'x_0', 'y_0', 'scale_x', 'scale_y'])
    temp = pd.merge(df, tile_origin, on='fov_id', how='left')

    cells_df = temp
    cells_df['x'] = temp.col + temp.x_0 * temp.scale_x  # col is the local x-coord. This line transforms local to global
    cells_df['y'] = temp.row + temp.y_0 * temp.scale_y  # row is the local y-coord. This line transforms local to global
    cells_df = cells_df.drop(columns=['x_0', 'y_0', 'scale_x', 'scale_y'])
    save_as = os.path.join(my_dir, 'all', 'expanded_cells.csv')
    cells_df.to_csv(save_as, index=False)

    logger.info('Done!')


def collate_label_image():
    ROOT_DIR = config.ROOT_DIR
    csv_root = os.path.join(ROOT_DIR, 'demo_data', 'human', 'label_image')
    target_file = os.path.join(csv_root, 'all', 'label_image.csv.gz')  # save as compressed gzip
    out = None
    tiles_across = 20
    tiles_down = 23
    tiles = []
    step = 0
    label_0 = 0
    for j in range(tiles_down):
        # loop over the rows and for each row stitch all csvs. The keep them in a list
        logger.info('Building row %d, / %d' % (j, tiles_down - 1))
        fovs = [int(i+j*tiles_across) for i in range(tiles_across)] # keep the fovs in this list
        csvs = [f"label_image_fov_{int(i+j*tiles_across):03d}.csv" for i in range(tiles_across)]
        filenames = [os.path.join(csv_root, csv) for csv in csvs]
        # df = pd.concat([pd.read_csv(f, index_col=False, header=None) for f in filenames], ignore_index=True, axis=1)
        df = [pd.read_csv(f, index_col=False, header=None) for f in filenames]
        max_labels = [d.max().max() for d in df]
        max_labels = accumulate(max_labels, lambda acc, elem: acc + elem if elem else 0)
        max_labels = list(max_labels)
        df_adj = [d[d > 0] + max_labels[i - 1] for i, d in enumerate(df)]  # adjust the labels
        df_adj = [d.fillna(0).astype(int) for d in df_adj]  # change the datatype
        df_adj = pd.concat(df_adj, ignore_index=True, axis=1)  # Now concatenate
        coo = coo_matrix(df_adj.values)
        coo.row += step
        tiles.append(coo)
        step = coo.shape[0] + step

    # find the correct dimensions for the final array
    dims = [x.shape for x in tiles]
    M = sum([d[0] for d in dims])
    N = list(set([d[1] for d in dims]))
    assert (len(N) == 1)
    N = N[0]
    # concatenate now the elements of the list
    _row = np.concatenate([x.row for x in tiles]).ravel().tolist()
    _col = np.concatenate([x.col for x in tiles]).ravel().tolist()
    _data = np.concatenate([x.data for x in tiles]).ravel().tolist()
    coo = coo_matrix((_data, (_row, _col)), shape=(M, N))
    out = pd.DataFrame({'row': _row, 'col': _col, 'label': _data})
    out.to_csv(target_file, index=False, compression='gzip')
    print(coo.data.size / (1024 ** 2))
    print(coo.shape)
    return coo


def adjust_label_images():
    '''
    reads each individual label_image and adjust the labels to a global context
    :return:
    '''
    ROOT_DIR = config.ROOT_DIR
    csv_root = os.path.join(ROOT_DIR, 'demo_data', 'human', 'label_image')
    filenames = os.listdir(csv_root)
    all_csv = [filename for filename in filenames if filename.endswith('csv') and filename.startswith('label_image_')]
    for fov in range(len(all_csv)):
        print('Doing fov %d' % fov)
        csv_str: str = f"label_image_fov_{int(fov):03d}.csv"
        csv_str = os.path.join(csv_root, csv_str)
        if fov == 0:
            cur = 0
        else:
            cur = df.max().max()

        df = pd.read_csv(csv_str, index_col=False, header=None)
        df = df[df != 0] + cur
        df = df.fillna(0).astype(int)
        df.to_csv(csv_str, header=False, index=False)

        print('ok')
    print('All done!')


def csv_splitter(target_file):
    '''
    Reads a single image_label file, splits it in to parts (every 1000 lines), then runs regionprops
    on the smaller files and then writes the expanded_cell flatilfes for each one. Then at then end
    it stiches together all the expanded_cells csv into one single flatfile
    :param target_file:
    :return:
    '''
    csvfile = open(target_file, 'r').readlines()
    _file = 0
    num_lines = len(csvfile)
    for i in range(num_lines):
        if i % 1000 == 0:
            filename = f"label_image_{int(_file):03d}.csv"
            filename = os.path.join('./', 'label_image', filename)
            open(str(filename), 'w+').writelines(csvfile[i:min(i+1000, num_lines)])
            _file += 1


def spot_in_cell():
    ROOT_DIR = os.path.join(config.ROOT_DIR, 'demo_data', 'human', 'spot_in_cell')
    spots_root = os.path.join(config.ROOT_DIR, 'demo_data', 'human', 'spots')
    labels_root = os.path.join(config.ROOT_DIR, 'demo_data', 'human', 'label_image')
    img_obj_root = os.path.join(config.ROOT_DIR, 'demo_data', 'human', 'expanded_cells')

    dir_files = os.listdir(spots_root)
    spots_csv = [filename for filename in dir_files if filename.endswith('.csv')]
    spots_fov = np.array([int(os.path.splitext(d)[0][-3:]) for d in spots_csv])

    dir_files = os.listdir(labels_root)
    labels_files = [filename for filename in dir_files if filename.endswith('.csv')]
    labels_fov = np.array([int(os.path.splitext(d)[0][-3:]) for d in labels_files])

    assert (np.all(spots_fov == labels_fov))

    _list = []
    for i in spots_fov:
        # i = 11
        logger.info('Doing fov %d / %d' % (i, len(spots_fov) - 1))
        res = pd.DataFrame(None, columns=['x_global',
                                          'y_global',
                                          'fov_id',
                                          'label_local',
                                          'x_cell',
                                          'y_cell'])
        missed = []
        at_bounds = []
        # save the output here
        target_file = os.path.join(ROOT_DIR, f'spot_in_cell_fov_{int(i):03d}.csv')

        # get the tile coords
        x_range, y_range, scaling_factor, tile_shape = get_tile_coords(i)

        lf = os.path.join(labels_root, f"label_image_fov_{int(i):03d}.csv")
        sf = os.path.join(spots_root, f"fov_{int(i):03d}.csv")
        cf = os.path.join(img_obj_root, f"expanded_cells_fov_{int(i):03d}.csv")

        spots_df = pd.read_csv(sf, index_col=False).dropna()

        if not spots_df.empty:

            mask = coofy(spots_df, [x_range[0], y_range[0]], scaling_factor, tile_shape)

            # attach an extra column with the global coordinates
            spots_df['x_global'] = spots_df.xc * scaling_factor
            spots_df['y_global'] = spots_df.yc * scaling_factor
            spots_df = spots_df \
                .astype({"x_global": int, "y_global": int}) \
                .sort_values(['x_global', 'y_global'], ascending=[True, True])

            # read the label_image as a dataframe
            my_labels = pd.read_csv(lf, index_col=False, header=None)  # O; background, 1: first cell, etc

            # logger.info('starting spot_map')
            sm = spot_map(mask, my_labels, tile_shape)
            # logger.info('spot_map ended')
            row = sm.row + y_range[0] * scaling_factor
            col = sm.col + x_range[0] * scaling_factor

            temp = pd.DataFrame({'x': col,
                             'y': row,
                             'label_local': sm.data,
                             'fov_id': i})

            cells_df = pd.read_csv(cf, index_col=False, header=0)
            cells_df['label_local'] = cells_df.cell_id + 1  # increase by one to align with the labels

            # attach an extra column with the global coordinates
            cells_df['x_cell'] = cells_df.col + x_range[0] * scaling_factor
            cells_df['y_cell'] = cells_df.row + y_range[0] * scaling_factor
            cells_df = cells_df.astype({"label_local": int, "fov_id": int})

            # attach now the correcponding cell centroid coordinates
            temp = temp.merge(cells_df, on=['label_local', 'fov_id'], how='left')

            res = pd.DataFrame({'x_global': temp.x,
                                'y_global': temp.y,
                                'fov_id': temp.fov_id,
                                'label_local': temp.label_local,
                                'x_cell': temp.x_cell,
                                'y_cell': temp.y_cell}).sort_values(['x_global', 'y_global'], ascending=[True, True])

            res = res.merge(pd.DataFrame({'x_global': spots_df.x_global,
                                          'y_global': spots_df.y_global,
                                          'target': spots_df.target}),
                            on=['x_global', 'y_global'],
                            how='left')

            # missed = [(x, y) for (x, y) in zip(spots_df.x_global.values, spots_df.y_global.values) if
            #  (x, y) not in zip(res.x_global.values, res.y_global.values)]
            #
            # at_bounds = [(x, y) for (x, y) in zip(spots_df.x_global.values, spots_df.y_global.values) if
            #  np.any([x == x_range[1] * scaling_factor, y == y_range[1] * scaling_factor])]

        # assert(sorted(missed) == sorted(at_bounds))

        # logger.info('Found %s spots of which %d line on the some' % (res.shape[0], res.dropna().shape[0]))
        _list.append(res)

        # write now to csv
        res.to_csv(target_file, index=False)

    # combine now everything and save the cvs
    df = pd.concat(_list, sort=False)

    df = df.sort_values(['x_global', 'y_global'], ascending=[True, True])
    logger.info('ok, start now')
    df = global_label(df)
    logger.info('done')

    df.to_csv(os.path.join(ROOT_DIR, 'all', 'spot_in_cell.csv.gz'), index=False, compression='gzip')


def global_label(df):
    centroids = sorted(set([x for x in zip(df.dropna().x_cell, df.dropna().y_cell)]))
    df['label_global'] = [centroids.index(d) + 1 if d[1] == d[1] else 0 for d in
                          zip(df.x_cell.values, df.y_cell.values)]
    return df


def global_label_2(df):
    # No idea why this one takes ages (20mins) to finish. It was fine yesterday.
    # It should take less than a sec
    # Even when run from a notebook works fine
    df['label_global'] = df.sort_values(['x_cell', 'y_cell'], na_position='first')[['x_cell', 'y_cell']]\
                             .fillna(0)\
                             .diff()\
                             .ne([0, 0])\
                             .any(1)\
                             .cumsum()-1
    return df


def spot_map(is_spot, image_label, dim):
    '''
    maps a spot to a label. If a spot has label==nan then it means the spot is on the background only
    If label = 1 then the spot lies within the first cell etc...
    :param is_spot:
    :param image_label:
    :param dim:
    :return:
    '''
    labeled_spots = coo_matrix(is_spot.toarray() * image_label)
    zero_list = [(i, j) for (i, j) in zip(is_spot.row, is_spot.col)
                 if (i, j) not in zip(labeled_spots.row, labeled_spots.col)]
    zero_row, zero_col = list(zip(*zero_list))
    zero_nan = np.nan * np.ones(len(zero_row))
    unlabeled_spots = coo_matrix((zero_nan, (zero_row, zero_col)), shape=dim+1)
    unlabeled_spots_arr = unlabeled_spots.toarray()
    unlabeled_spots = coo_matrix(unlabeled_spots_arr[:-1, :-1])
    row = np.append(labeled_spots.row, unlabeled_spots.row)
    col = np.append(labeled_spots.col, unlabeled_spots.col)
    data = np.append(labeled_spots.data, unlabeled_spots.data)
    assert(np.all(np.sort(row) == np.sort(is_spot.row)))
    assert(np.all(np.sort(col) == np.sort(is_spot.col)))
    out = coo_matrix((data, (row, col)), shape=image_label.shape)
    return out


def coofy(df, origin, scaling_factor, tile_shape):
    '''
    makes all coords local and returns a sparse array
    :param df:
    :param origin:
    :param scaling_factor:
    :param tile_shape:
    :return:
    '''

    x0 = origin[0]
    y0 = origin[1]
    x = (df.xc - x0) * scaling_factor
    y = (df.yc - y0) * scaling_factor

    ones = np.ones(x.shape)
    coo = coo_matrix((ones, (y, x)), shape=tile_shape + 1)  # Why the +1? starfish includes the boundaries,
    # hence the extra pixel. That sounds like a bug maybe? a tile of size 2000-by-2000 should have a range
    # [0, 1999] going across and [0, 1999] going down
    coo_arr = coo.toarray()
    out = coo_matrix(coo_arr[:-1, :-1])  # why the -1? See comment right above
    return out


def set_labels(df):
    '''
    assigns the labels
    :param df:
    :return:
    '''
    my_set = set(list(zip(df.x_cell, df.y_cell)))
    ordered_set = sorted(my_set, key=lambda tup: (tup[0], tup[1]))
    cell_df = pd.DataFrame(ordered_set, columns=['x_cell', 'y_cell'])\
        .reset_index()\
        .rename(columns={'index': 'label'})

    # I want my labels to start from 1
    cell_df['label'] = cell_df.label + 1

    # for each spot stick now the coords and the label of the closest cell. That applies only to
    # spots whose coordinates overlap with a label from the relavant label_image.
    out = df.merge(cell_df, on=['x_cell', 'y_cell'], how='left')

    # cell_id is local id wrt to the relevant fov. No need to keep it. Drop it
    out = out.drop(['cell_id'], axis=1)\
        .sort_values(['x_cell', 'y_cell'], ascending=[True, True])
    return out


def get_tile_coords(fov):
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
    base64string = base64.b64encode(bytes('%s:%s' % (credentials.USER, credentials.PASSWD), 'ascii'))
    request.add_header("Authorization", "Basic %s" % base64string.decode('utf-8'))
    with urllib.request.urlopen(request) as result:
        data = json.loads(result.read().decode())
        # data = json.load(result) # i think that also works, dont know whats best!
        xc = data['tiles'][0]['coordinates']['xc']
        yc = data['tiles'][0]['coordinates']['yc']
        default_tile_shape = np.array(data['default_tile_shape'])
        tile_shape = np.array(data['tiles'][0]['tile_shape'])
        scale_x = data['default_tile_shape'][0] / np.diff(data['tiles'][0]['coordinates']['xc'])[0]
        scale_y = data['default_tile_shape'][1] / np.diff(data['tiles'][0]['coordinates']['yc'])[0]
        assert(scale_y == scale_x)
        assert(np.all(default_tile_shape == tile_shape))
        scaling_factor = scale_x

    return xc, yc, scaling_factor, tile_shape


def read_single_cell():
    '''
    reads expression table (columns: cell types, rows: gene names) and returns a dataframe
    of dimension numbers of genes-by-number of unique cell types describing the mean expression
    of each gene per cell type.
    Cell types are assumed to be in the first line of the csv and gene names in the first column
    :return:
    '''

    # NOTE: maybe set that path in the config file
    csv_path = os.path.join(config.ROOT_DIR, 'demo_data', 'human', 'scRNA', 'scRNAseqHCA120GenePanel.csv')
    scRNA_df = pd.read_csv(csv_path, index_col=0, header=0)

    # group by column names and find the meam
    out = scRNA_df.groupby(scRNA_df.columns, axis=1).mean()

    return out



