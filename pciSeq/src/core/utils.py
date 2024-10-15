import sys
import tempfile, shutil
from tqdm import tqdm
from urllib.parse import urlparse
from urllib.request import urlopen
import numpy as np
import pandas as pd
import pickle
import dask
import os
import glob
import subprocess
import numbers
from email.parser import BytesHeaderParser
from scipy.spatial import cKDTree
import webbrowser
import logging

utils_logger = logging.getLogger(__name__)


def negBinLoglik(x, r, p):
    '''
    Negative Binomial loglikehood
    :param x:
    :param r:
    :param p:
    :return:
    '''

    # sanity check
    # assert (np.all(da_x.coords['cell_id'].data == da_p.coords['cell_id'])), 'gene counts and beta probabilities are not aligned'
    # assert (np.all(da_x.coords['gene_name'].data == da_p.coords['gene_name'])), 'gene counts and beta probabilities are not aligned'

    nC, nG = x.shape
    x = x[:, :, None]
    nK = p.shape[-1]
    # contr = np.zeros([nC, nG, nK])
    contr = x * np.log(p) + r * np.log(1 - p)
    # ne.evaluate("x * log(p) + r * log(1 - p)", out=contr)
    return contr


def softmax(X, theta=1.0, axis=None):
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
    y = y - np.expand_dims(np.max(y, axis=axis), axis)

    # exponentiate y
    y = np.exp(y)

    # take the sum along the specified axis
    ax_sum = np.expand_dims(np.sum(y, axis=axis), axis)

    # finally: divide elementwise
    p = y / ax_sum

    # flatten if X was 1D
    if len(X.shape) == 1: p = p.flatten()

    return p


def hasConverged(spots, p0, tol):
    p1 = spots.parent_cell_prob
    if p0 is None:
        p0 = np.zeros(p1.shape)
    delta = np.max(np.abs(p1 - p0))
    converged = (delta < tol)
    return converged, delta


def splitter_mb(filepath, mb_size):
    """ Splits a text file in (almost) equally sized parts on the disk. Assumes that there is a header in the first line
    :param filepath: The path of the text file to be broken up into smaller files
    :param mb_size: size in MB of each chunk
    :return:
    """
    handle = open(filepath, 'r')
    OUT_DIR = os.path.join(os.path.splitext(filepath)[0] + '_split')

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    else:
        files = glob.glob(OUT_DIR + '/*.*')
        for f in files:
            os.remove(f)

    n = 0
    size = None
    header_line = next(handle)
    file_out, handle_out = _get_file(OUT_DIR, filepath, n, header_line)
    for line in handle:
        size = os.stat(file_out).st_size
        if size > mb_size * 1024 * 1024:
            utils_logger.info('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
            n += 1
            handle_out.close()
            file_out, handle_out = _get_file(OUT_DIR, filepath, n, header_line)
        handle_out.write(str(line))

    # print(str(file_out) + " file size = \t" + str(size))
    print('saved %s with file size %4.3f MB' % (file_out, size / (1024 * 1024)))
    handle_out.close()


def splitter_n(filepath, n):
    """ Splits a text file in n smaller files
    :param filepath: The path of the text file to be broken up into smaller files
    :param n: determines how many smaller files will be created
    :return:
    """
    filename_ext = os.path.basename(filepath)
    [filename, ext] = filename_ext.split('.')

    OUT_DIR = os.path.join(os.path.splitext(filepath)[0] + '_split')

    if ext == 'json':
        df = pd.read_json(filepath)
    elif ext == 'tsv':
        df = pd.read_csv(filepath, sep='\t')
    else:
        df = None

    df_list = np.array_split(df, n)
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    else:
        files = glob.glob(OUT_DIR + '/*.' + ext)
        for f in files:
            os.remove(f)

    for i, d in enumerate(df_list):
        fname = os.path.join(OUT_DIR, filename + '_%d.%s' % (i, ext))
        if ext == 'json':
            d.to_json(fname, orient='records')
        elif ext == 'tsv':
            d.to_csv(fname, sep='\t', index=False)


def _get_file(OUT_DIR, filepath, n, header_line):
    [filename, ext] = os.path.basename(filepath).split('.')
    file = os.path.join(OUT_DIR, filename + '_%d.%s' % (n, ext))
    handle = open(file, "a")
    handle.write(header_line)
    return file, handle


def download_url_to_file(url, dst, progress=True):
    r"""Download object at the given URL to a local path.
            Thanks to torch, slightly modified
    Args:
        url (string): URL of the object to download
        dst (string): Full path where object will be saved, e.g. `/tmp/temporary_file`
        progress (bool, optional): whether or not to display a progress bar to stderr
            Default: True
    """
    file_size = None
    u = urlopen(url)
    meta = u.info()
    if hasattr(meta, 'getheaders'):
        content_length = meta.getheaders("Content-Length")
    else:
        content_length = meta.get_all("Content-Length")
    if content_length is not None and len(content_length) > 0:
        file_size = int(content_length[0])
    # We deliberately save it in a temp file and move it after
    dst = os.path.expanduser(dst)
    dst_dir = os.path.dirname(dst)
    f = tempfile.NamedTemporaryFile(delete=False, dir=dst_dir)
    try:
        with tqdm(total=file_size, disable=not progress,
                  unit='B', unit_scale=True, unit_divisor=1024) as pbar:
            while True:
                buffer = u.read(8192)
                if len(buffer) == 0:
                    break
                f.write(buffer)
                pbar.update(len(buffer))
        f.close()
        shutil.move(f.name, dst)
    finally:
        f.close()
        if os.path.exists(f.name):
            os.remove(f.name)


def load_from_url(url):
    # files = []
    parts = urlparse(url)
    filename = os.path.basename(parts.path)
    if not os.path.exists(filename):
        sys.stderr.write('Downloading: "{}" to {}\n'.format(url, filename))
        download_url_to_file(url, filename)
    return filename


def get_out_dir(path=None, sub_folder=''):
    if path is None or path[0] == 'default':
        out_dir = os.path.join(tempfile.gettempdir(), 'pciSeq')
    else:
        out_dir = os.path.join(path[0], sub_folder, 'pciSeq')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    return out_dir


def get_pciSeq_install_dir():
    p = subprocess.run(['pip', 'show', 'pciSeq'], stdout=subprocess.PIPE)
    h = BytesHeaderParser().parsebytes(p.stdout)
    assert h['Location'] is not None, 'Could not locate pciSeq installation folder, maybe the package is not installed.'
    return os.path.join(h['Location'], 'pciSeq')


def gaussian_contour(mu, cov, sdwidth=3):
    """
    Draws an ellipsoid for a given covariance matrix cov
    and mean vector mu

    Example
    cov_1 = [[1, 0.5], [0.5, 1]]
    means_1 = [1, 1]

    cov_2 = [[1, -0.7], [-0.7, 1]]
    means_2 = [2, 1.5]

    cov_3 = [[1, 0], [0, 1]]
    means_3 = [0, 0]

    ellipsis_1 = gaussian_contour(means_1, cov_1)
    ellipsis_2 = gaussian_contour(means_2, cov_2)
    ellipsis_3 = gaussian_contour(means_3, cov_3)
    ellipsis_3b = gaussian_contour(means_3, cov_3, sdwidth=2)
    ellipsis_3c = gaussian_contour(means_3, cov_3, sdwidth=3)

    plt.plot(ellipsis_1[0], ellipsis_1[1], c='b')
    plt.plot(ellipsis_2[0], ellipsis_2[1], c='r')
    plt.plot(ellipsis_3[0], ellipsis_3[1], c='g')
    plt.plot(ellipsis_3b[0], ellipsis_3b[1], c='g')
    plt.plot(ellipsis_3c[0], ellipsis_3c[1], c='g')
    """

    # cov_00 = sigma_x * sigma_x
    # cov_10 = rho * sigma_x * sigma_y
    # cov_11 = sigma_y * sigma_y
    # cov = np.array([[cov_00, cov_10], [cov_10, cov_11]])
    mu = np.array(mu)

    npts = 40
    tt = np.linspace(0, 2 * np.pi, npts)
    ap = np.zeros((2, npts))
    x = np.cos(tt)
    y = np.sin(tt)
    ap[0, :] = x
    ap[1, :] = y

    eigvals, eigvecs = np.linalg.eig(cov)
    eigvals = sdwidth * np.sqrt(eigvals)
    eigvals = eigvals * np.eye(2)

    vd = eigvecs.dot(eigvals)
    out = vd.dot(ap) + mu.reshape(2, -1)

    return np.array(list(zip(*out)), dtype=np.float32)


def read_image_objects(img_obj, cfg):
    meanCellRadius = np.mean(np.sqrt(img_obj.area / np.pi)) * 0.5
    relCellRadius = np.sqrt(img_obj.area / np.pi) / meanCellRadius

    # append 1 for the misreads
    relCellRadius = np.append(1, relCellRadius)

    InsideCellBonus = cfg['InsideCellBonus']
    if not InsideCellBonus:
        # This is more for clarity. The operation below will work fine even if InsideCellBonus is False
        InsideCellBonus = 0

    # if InsideCellBonus == 0 then CellAreaFactor will be equal to 1.0
    numer = np.exp(-relCellRadius ** 2 / 2) * (1 - np.exp(InsideCellBonus)) + np.exp(InsideCellBonus)
    denom = np.exp(-0.5) * (1 - np.exp(InsideCellBonus)) + np.exp(InsideCellBonus)
    CellAreaFactor = numer / denom

    out = {}
    out['area_factor'] = CellAreaFactor.astype(np.float32)
    out['rel_radius'] = relCellRadius.astype(np.float32)
    out['area'] = np.append(np.nan, img_obj.area.astype(np.uint32))
    out['x0'] = np.append(-sys.maxsize, img_obj.x0.values).astype(np.float32)
    out['y0'] = np.append(-sys.maxsize, img_obj.y0.values).astype(np.float32)
    out['z0'] = np.append(-sys.maxsize, img_obj.z0.values).astype(np.float32)
    out['cell_label'] = np.append(0, img_obj.label.values).astype(np.uint32)
    if 'old_label' in img_obj.columns:
        out['cell_label_old'] = np.append(0, img_obj.old_label.values).astype(np.uint32)
    # First cell is a dummy cell, a super neighbour (ie always a neighbour to any given cell)
    # and will be used to get all the misreads. It was given the label=0 and some very small
    # negative coords

    return out, meanCellRadius.astype(np.float32)


def keep_labels_unique(scdata):
    """
    In the single cell data you might find cases where two or more rows have the same gene label
    In these cases keep the row with the highest total gene count
    """

    # 1. get the row total and assign it to a new column
    scdata = scdata.assign(total=scdata.sum(axis=1))

    # 2. rank by gene label and total gene count and keep the one with the highest total
    scdata = scdata.sort_values(['gene_name', 'total'], ascending=[True, False]).groupby('gene_name').head(1)

    # 3. Drop the total column and return
    return scdata.drop(['total'], axis=1)


# @dask.delayed
def scaled_exp(cell_area_factor, sc_mean_expressions, inefficiency):
    if np.all(cell_area_factor == 1):
        subscripts = 'gk, g -> gk'
        operands = [sc_mean_expressions, inefficiency]
    else:
        subscripts = 'c, gk, g -> cgk'
        operands = [cell_area_factor, sc_mean_expressions, inefficiency]
    out = np.einsum(subscripts, *operands)
    return out


def adjust_for_anisotropy(spots, voxel_size):
    gene_col = spots.gene_name.values[:, None]
    z_plane = np.float32(spots.z_plane.values[:, None])

    data = spots[['x', 'y', 'z_plane']]
    spots_adj = anisotropy_calc(data, voxel_size)

    spots_adj = np.hstack([gene_col, spots_adj, z_plane])
    spots_adj = pd.DataFrame(spots_adj, columns=['gene_name', 'x', 'y', 'z', 'z_plane'])
    spots_adj = spots_adj.astype({'x': 'float32',
                                  'y': 'float32',
                                  'z': 'float32',
                                  'z_plane': 'float32',
                                  'gene_name': str})
    return spots_adj  # Nspots x 3


def anisotropy_calc(data, voxel_size):
    x, y, z = voxel_size
    Sx = x / x
    Sy = y / x
    Sz = z / x
    scaling_matrix = np.array([
        [Sx, 0, 0],
        [0, Sy, 0],
        [0, 0, Sz]
    ], dtype=np.float32)
    out = scaling_matrix.dot(data.T)  # 3 x Nspots
    return out.T  # Nspots x 3


def truncate_data(spots, label_image, z_min, z_max):
    if z_min is not None and z_max is not None:
        utils_logger.info(' Truncating masks and spots. Keeping those between frame %d and %d only' % (z_min, z_max))
        spots = truncate_spots(spots, z_min, z_max)
    return spots, label_image


def truncate_spots(spots, zmin, zmax):
    spots_min = spots[(spots.z_stack <= zmax) & (spots.z_stack >= zmin)]
    spots_min = spots_min.assign(z_stack=spots_min.z_stack - zmin)
    # out = spots_min.z - zmin
    return spots_min.reset_index(drop=True)


def get_img_shape(coo):
    n = len(coo)
    img_shape = set([d.shape for d in coo])
    assert len(img_shape) == 1, 'pages do not have the same shape'
    img_shape = img_shape.pop()
    w = img_shape[1]
    h = img_shape[0]
    return [n, h, w]


def gaussian_ellipsoid_props(cov, sdwidth=3):
    """
    get the scaling, rotation of the ellipsoid
    """
    tol = 1.0e-10
    cov = np.where(cov < tol, 0, cov)
    eigvals, eigvecs = np.linalg.eig(cov)
    scaling = sdwidth * np.sqrt(eigvals)
    # rotation = roll_pitch_yaw(eigvecs)
    rotation = euler_angles(eigvecs.T)
    scaling = scaling.tolist()
    return scaling, rotation


def euler_angles(r):
    theta_x = np.arctan2(r[2, 1], r[2, 2])
    theta_y = np.arcsin(r[2, 0])
    theta_z = np.arctan2(r[1, 0], r[0, 0])

    return [theta_x, theta_y, theta_z]


def gaussian_contour(mu, cov, sdwidth=3):
    """
    Draws an ellipsoid for a given covariance matrix cov
    and mean vector mu

    Example
    cov_1 = [[1, 0.5], [0.5, 1]]
    means_1 = [1, 1]

    cov_2 = [[1, -0.7], [-0.7, 1]]
    means_2 = [2, 1.5]

    cov_3 = [[1, 0], [0, 1]]
    means_3 = [0, 0]

    ellipsis_1 = gaussian_contour(means_1, cov_1)
    ellipsis_2 = gaussian_contour(means_2, cov_2)
    ellipsis_3 = gaussian_contour(means_3, cov_3)
    ellipsis_3b = gaussian_contour(means_3, cov_3, sdwidth=2)
    ellipsis_3c = gaussian_contour(means_3, cov_3, sdwidth=3)

    plt.plot(ellipsis_1[0], ellipsis_1[1], c='b')
    plt.plot(ellipsis_2[0], ellipsis_2[1], c='r')
    plt.plot(ellipsis_3[0], ellipsis_3[1], c='g')
    plt.plot(ellipsis_3b[0], ellipsis_3b[1], c='g')
    plt.plot(ellipsis_3c[0], ellipsis_3c[1], c='g')
    """

    # cov_00 = sigma_x * sigma_x
    # cov_10 = rho * sigma_x * sigma_y
    # cov_11 = sigma_y * sigma_y
    # cov = np.array([[cov_00, cov_10], [cov_10, cov_11]])
    mu = np.array(mu)

    npts = 40
    tt = np.linspace(0, 2 * np.pi, npts)
    ap = np.zeros((2, npts))
    x = np.cos(tt)
    y = np.sin(tt)
    ap[0, :] = x
    ap[1, :] = y

    eigvals, eigvecs = np.linalg.eig(cov)
    eigvals = sdwidth * np.sqrt(eigvals)
    eigvals = eigvals * np.eye(2)

    vd = eigvecs.dot(eigvals)
    out = vd.dot(ap) + mu.reshape(2, -1)

    return np.array(list(zip(*out)), dtype=np.float32)


def recover_original_labels(cellData, geneData, remapping):
    labels_dict = dict(zip(remapping.values(), remapping.keys()))
    cellData = cellData.assign(Cell_Num=cellData.Cell_Num.map(lambda x: labels_dict[x]))

    # geneData = geneData.assign(neighbour=geneData.neighbour.map(lambda x: swaped_dict[x]))
    geneData = geneData.assign(neighbour=geneData.neighbour.map(lambda x: fetch_label(x, labels_dict)))
    geneData = geneData.assign(neighbour_array=geneData.neighbour_array.map(lambda x: fetch_label(x, labels_dict)))
    return cellData, geneData


def fetch_label(x, d):
    x = [x] if isinstance(x, numbers.Number) else x
    out = [d[v] for v in x]
    return out[0] if len(out) == 1 else out


def serialise(varBayes, debug_dir):
    if not os.path.exists(debug_dir):
        os.makedirs(debug_dir)
    pickle_dst = os.path.join(debug_dir, 'pciSeq.pickle')
    with open(pickle_dst, 'wb') as outf:
        pickle.dump(varBayes, outf)
        utils_logger.info('Saved at %s' % pickle_dst)


def purge_spots(spots, sc):
    drop_spots = list(set(spots.gene_name) - set(sc.index))
    utils_logger.warning('Found %d genes that are not included in the single cell data' % len(drop_spots))
    idx = ~ np.in1d(spots.gene_name, drop_spots)
    spots = spots.iloc[idx]
    utils_logger.warning('Removed from spots: %s' % drop_spots)
    return spots


def export_db_tables(out_dir, con):
    tables = con.get_db_tables()
    for table in tables:
        export_db_table(table, out_dir, con)


def export_db_table(table_name, out_dir, con):
    df = con.from_redis(table_name)
    fname = os.path.join(out_dir, table_name + '.csv')
    df.to_csv(fname, index=False)
    utils_logger.info('Saved at %s' % fname)


def get_closest(spots, query_vals, voxel_size=None):

    # Check if the input is a DataFrame
    if isinstance(query_vals, pd.DataFrame):
        pass

    # If input is a list, check if it's a list of dictionaries
    elif isinstance(query_vals, list) and all(isinstance(item, dict) for item in query_vals):
        query_vals = pd.DataFrame(query_vals)  # Input is already a list of dicts, return as is

    if voxel_size is None:
        voxel_size = [1, 1, 1]

    # Filter spots based on gene_name
    # query_vals = pd.DataFrame(query_vals)
    if 'z_plane' in query_vals.columns:
        query_vals = adjust_for_anisotropy(query_vals, voxel_size)

    gene_names = query_vals['gene_name'].tolist()

    mask = np.isin(gene_names, spots['gene_name'].values)
    if not np.all(mask):
        missing_genes = np.array(gene_names)[~mask]  # Get the missing genes
        # Format and print the message
        print(f"Genes {', '.join(missing_genes)} are not in the neighborhood. Skipping.")

    mask = spots['gene_name'].isin(gene_names)
    filtered_spots = spots[mask]
    if filtered_spots.empty:
        result = None
    else:
        # Create KDTree from filtered spots
        tree = cKDTree(filtered_spots[['x', 'y', 'z']])

        # Find nearest neighbors
        distances, indices = tree.query(np.column_stack((query_vals['x'], query_vals['y'], query_vals['z'])))

        # Create result DataFrame
        result = filtered_spots.iloc[indices].reset_index()
        result = result[result['gene_name'] == query_vals['gene_name']]

        # Sort by spot_id (index of the original spots DataFrame)
        result = result.sort_values('spot_id')

        # Select and rename columns
        result = result[['gene_name', 'spot_id', 'x', 'y', 'z']]

    return result


def gene_loglik_contributions_scatter(data, filename='interactive_scatter.html'):
    # Convert data and classes to JSON
    # data_json = json.dumps(data)
    # classes_json = json.dumps(classes)

    # HTML and JavaScript code
    html_code = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Interactive Scatter Plot</title>
        <script src="https://d3js.org/d3.v7.min.js"></script>
        <style>
            body {{
                margin: 0;
                padding: 0;
                height: 100vh;
                display: flex;
                flex-direction: column;
            }}
            #top-space {{
               height: 40px;
               display: flex;
               justify-content: space-between;
               align-items: center;
               padding: 0 10px;
            }}
            #plot-area {{
                flex-grow: 1;
                width: 100%;
                position: relative;
            }}
            .tooltip {{
                position: absolute;
                background-color: rgba(0, 0, 0, 0.7);
                color: white;
                padding: 5px;
                border-radius: 5px;
                font-size: 12px;
                pointer-events: none;
            }}
            .axis-label {{
                font-size: 14px;
                font-weight: bold;
            }}
            .plot-title {{
                font-size: 18px;
                font-weight: bold;
                text-anchor: middle;
            }}
            .plot-subtitle {{
                font-size: 14px;
                text-anchor: middle;
            }}
            #dropdown {{
                font-size: 14px;
                padding: 5px;
            }}
        </style>
    </head>
    <body>
        <div id="top-space">
            <!-- <select id="dropdown"></select> -->
            <div id="title-space"></div>
        </div>
        <div id="plot-area"></div>
        <script>
            const data = {data['contr']};
            const classes = {data['class_names']};
            const gene_names = {data['gene_names']}
            let currentUserClass = "{data['user_class']}";
            let currentAssignedClass = "{data['assigned_class']}";

            // Populate dropdown
            const dropdown = d3.select('#dropdown');
            dropdown.selectAll('option')
                .data(classes)
                .enter()
                .append('option')
                .text(d => d)
                .attr('value', d => d);

            const margin = {{top: 60, right: 80, bottom: 50, left: 100}};
            const width = window.innerWidth - margin.left - margin.right;
            const height = window.innerHeight * 0.34 - margin.top - margin.bottom;

            const svg = d3.select('#plot-area')
                .append('svg')
                .attr('width', width + margin.left + margin.right)
                .attr('height', height + margin.top + margin.bottom)
                .append('g')
                .attr('transform', `translate(${{margin.left}},${{margin.top}})`);

            const x = d3.scaleLinear()
                .domain(d3.extent(data[currentAssignedClass]))
                .range([0, width]);

            const y = d3.scaleLinear()
                .domain(d3.extent(data[currentUserClass]))
                .range([height, 0]);

            const xAxis = svg.append('g')
                .attr('class', 'x-axis')
                .attr('transform', `translate(0,${{height}})`)
                .call(d3.axisBottom(x));

            const yAxis = svg.append('g')
                .attr('class', 'y-axis')
                .call(d3.axisLeft(y));

            // Add plot title
            svg.append("text")
                .attr("class", "plot-title")
                .attr("x", width / 2)
                .attr("y", -margin.top / 2)
                .text("Loglikelihood contributions for cell: {data['cell_num']}");

            // Add plot subtitle
            svg.append("text")
                .attr("class", "plot-subtitle")
                .attr("x", width / 2)
                .attr("y", -margin.top / 2 + 20)
                .text("Assigned class: {data['assigned_class']}");

            // Add x-axis label
            svg.append("text")
                .attr("class", "axis-label x-axis-label")
                .attr("x", width / 2)
                .attr("y", height + margin.bottom - 10)
                .style("text-anchor", "middle")
                .text("{data['assigned_class']}");

            // Add y-axis label
            svg.append("text")
                .attr("class", "axis-label y-axis-label")
                .attr("transform", "rotate(-90)")
                .attr("x", -height / 2)
                .attr("y", -margin.left + 40)
                .style("text-anchor", "middle")
                .text("{data['user_class']}");

            // Create tooltip
            const tooltip = d3.select("body").append("div")
                .attr("class", "tooltip")
                .style("opacity", 0);

            let myData = data[currentAssignedClass].map((x, index) => ({{
                'x': x,
                'y': data[currentUserClass][index],
                'name': gene_names[index]
            }}));

            function doData(data, currentAssignedClass, currentUserClass){{
                return data[currentAssignedClass].map((x, index) => ({{
                    'x': x,
                    'y': data[currentUserClass][index]
                }}));
            }}

            const dots = svg.selectAll('circle')
                .data(myData)
                .enter()
                .append('circle')
                .attr('cx', d => x(d.x))
                .attr('cy', d => y(d.y))
                .attr('r', 5)
                .style('fill', 'blue')
                .on("mouseover", function(event, d) {{
                    tooltip.transition()
                        .duration(200)
                        .style("opacity", .9);
                    tooltip.html(`Name: ${{d.name}}<br>X: ${{d.x.toFixed(3)}}<br>Y: ${{d.y.toFixed(3)}}`)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                }})
                .on("mousemove", function(event, d) {{
                    tooltip.style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                }})
                .on("mouseout", function(d) {{
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
                }});


            // Add interpretation guide
            const guideText = [
                "Interpretation Guide:",
                `• Genes on diagonal: Contribute equally to both cell types`,
                `• Genes above diagonal: Support classification as {data['user_class']}`,
                `• Genes below diagonal: Support classification as {data['assigned_class']}`,
                `• Distance from diagonal: Strength of support for one type over the other`
            ];

            const guide = svg.append("g")
                .attr("class", "interpretation-guide")
                .attr("transform", `translate(${{width - 10}}, ${{height - 10}})`);

            guide.selectAll("text")
                .data(guideText)
                .enter()
                .append("text")
                .attr("x", 0)
                .attr("y", (d, i) => i * 15)
                .style("text-anchor", "start")
                .style("font-size", "12px")
                .text(d => d);

            // Add a semi-transparent background to the guide
            const guideBBox = guide.node().getBBox();
            guide.insert("rect", ":first-child")
                .attr("x", guideBBox.x - 5)
                .attr("y", guideBBox.y - 5)
                .attr("width", guideBBox.width + 10)
                .attr("height", guideBBox.height + 10)
                .attr("fill", "rgba(255, 223, 186, 0.7)");  // Light pastel orange


            // Adjust the position of the guide to the bottom right
            guide.attr("transform", `translate(${{width - guideBBox.width - 15}}, ${{height - guideBBox.height - 15}})`);


            // Add diagonal line (y=x)
            const diagonalLine = svg.append('line')
                .attr('class', 'diagonal-line')
                .attr('x1', x(Math.min(x.domain()[0], y.domain()[0])))
                .attr('y1', y(Math.min(x.domain()[0], y.domain()[0])))
                .attr('x2', x(Math.max(x.domain()[1], y.domain()[1])))
                .attr('y2', y(Math.max(x.domain()[1], y.domain()[1])))
                .attr('stroke', 'red')
                .attr('stroke-width', 2)
                .attr('stroke-dasharray', '5,5');

            // Add zoom and pan handling
            let zoom = 1;
            let translateX = 0;
            let translateY = 0;
            let isDragging = false;
            let startX, startY;

            function updatePlot() {{
                const newX = x.copy().range([translateX, width * zoom + translateX]);
                const newY = y.copy().range([height * zoom + translateY, translateY]);

                // Update the axes
                xAxis.call(d3.axisBottom(newX));
                yAxis.call(d3.axisLeft(newY));

                // Update the position of the dots
                dots.attr('cx', d => newX(d.x))
                    .attr('cy', d => newY(d.y));

                // Update the diagonal line
                const minDomain = Math.min(newX.domain()[0], newY.domain()[0]);
                const maxDomain = Math.max(newX.domain()[1], newY.domain()[1]);
                diagonalLine
                    .attr('x1', newX(minDomain))
                    .attr('y1', newY(minDomain))
                    .attr('x2', newX(maxDomain))
                    .attr('y2', newY(maxDomain));
            }}

            d3.select('#plot-area').on('wheel', function(event) {{
                event.preventDefault();
                const delta = event.deltaY;
                if (delta > 0) {{
                    zoom *= 0.9;
                }} else {{
                    zoom *= 1.1;
                }}
                zoom = Math.max(0.1, Math.min(zoom, 10));
                updatePlot();
            }});

            d3.select('#plot-area')
                .on('mousedown', function(event) {{
                    isDragging = true;
                    startX = event.clientX - translateX;
                    startY = event.clientY - translateY;
                }})
                .on('mousemove', function(event) {{
                    if (isDragging) {{
                        translateX = event.clientX - startX;
                        translateY = event.clientY - startY;
                        updatePlot();
                    }}
                }})
                .on('mouseup', function() {{
                    isDragging = false;
                }})
                .on('mouseleave', function() {{
                    isDragging = false;
                }});

            // Resize function
            function resizePlot() {{
                const newWidth = window.innerWidth - margin.left - margin.right;
                const newHeight = window.innerHeight * 0.34 - margin.top - margin.bottom;

                svg.attr('width', newWidth + margin.left + margin.right)
                   .attr('height', newHeight + margin.top + margin.bottom);

                x.range([0, newWidth]);
                y.range([newHeight, 0]);

                xAxis.attr('transform', `translate(0,${{newHeight}})`).call(d3.axisBottom(x));
                yAxis.call(d3.axisLeft(y));

                dots.attr('cx', d => x(d.x))
                    .attr('cy', d => y(d.y));

                const minDomain = Math.min(x.domain()[0], y.domain()[0]);
                const maxDomain = Math.max(x.domain()[1], y.domain()[1]);
                diagonalLine
                    .attr('x1', x(minDomain))
                    .attr('y1', y(minDomain))
                    .attr('x2', x(maxDomain))
                    .attr('y2', y(maxDomain));

                // Update axis labels, title, and subtitle
                svg.select(".plot-title")
                    .attr("x", newWidth / 2);
                svg.select(".plot-subtitle")
                    .attr("x", newWidth / 2);
                svg.select(".x-axis-label")
                    .attr("x", newWidth / 2)
                    .attr("y", newHeight + margin.bottom - 10);
                svg.select(".y-axis-label")
                    .attr("x", -newHeight / 2);
            }}

            // Add event listener for window resize
            window.addEventListener('resize', resizePlot);

            // Update function for changing the plot
            function updatePlotForNewClass(newUserClass) {{
                currentUserClass = newUserClass;

                var data_upd = doData(data, currentAssignedClass, currentUserClass)

                // Update y-axis label
                svg.select(".y-axis-label")
                    .text(currentUserClass);

                // Update plot subtitle
                svg.select(".plot-subtitle")
                    .text(`Assigned class: ${{currentAssignedClass}} vs Selected class: ${{currentUserClass}}`);

                // Update interpretation guide
                guide.selectAll("text")
                    .data([
                        "Interpretation Guide:",
                        `• Genes on diagonal: Contribute equally to both cell types`,
                        `• Genes above diagonal: Support classification as ${{currentUserClass}}`,
                        `• Genes below diagonal: Support classification as ${{currentAssignedClass}}`,
                        `• Distance from diagonal: Strength of support for one type over the other`
                    ])
                    .text(d => d);

                // Update y-axis scale
                const newYExtent = d3.extent(data[currentUserClass]);
                y.domain(newYExtent);
                yAxis.transition().duration(1000).call(d3.axisLeft(y));

                // Update dots
                dots.transition()
                    .duration(1000)
                    .attr('cy', d => 1.0)

                // Update diagonal line
                const minDomain = Math.min(x.domain()[0], y.domain()[0]);
                const maxDomain = Math.max(x.domain()[1], y.domain()[1]);
                diagonalLine.transition()
                    .duration(1000)
                    .attr('y1', y(minDomain))
                    .attr('y2', y(maxDomain));

                // Update tooltip
                dots.on("mouseover", function(event, d) {{
                     tooltip.transition()
                         .duration(200)
                         .style("opacity", .9);
                     tooltip.html(`Name: ${{d.name}}<br>X: ${{d[currentAssignedClass].toFixed(3)}}<br>Y: ${{d[currentUserClass].toFixed(3)}}`)
                         .style("left", (event.pageX + 10) + "px")
                         .style("top", (event.pageY - 28) + "px");
                 }});
            }}

            // Dropdown event listener
            d3.select('#dropdown').on('change', function() {{
                const selectedValue = d3.select(this).property('value');
                console.log('Selected value:', selectedValue);
                updatePlotForNewClass(selectedValue);
            }});

        </script>
    </body>
    </html>
    """

    # Save the HTML to a file
    with open(filename, 'w') as f:
        f.write(html_code)

    # Open the HTML file in the default web browser
    webbrowser.open('file://' + os.path.realpath(filename))
