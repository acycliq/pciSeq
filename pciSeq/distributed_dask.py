import pandas as pd
import numpy  as np
import dask
import dask.array as da
import dask.dataframe as dd
from typing import Tuple
from scipy.sparse import load_npz, coo_matrix
from pciSeq import fit as fit_chunk
import time


import os


def filter_spots(spots, bbox):
    """
    removes out of bounds spots (if any...)
    """
    ymin = bbox[0][0]
    ymax = bbox[0][1]
    xmin = bbox[1][0]
    xmax = bbox[1][1]
    mask_x = (spots.x >= xmin) & (spots.x < xmax)
    mask_y = (spots.y >= ymin) & (spots.y < ymax)
    df = spots[mask_x & mask_y]
    df = df.assign(x=df.x-xmin)
    df = df.assign(y=df.y-ymin)
    origin = (xmin, ymin)
    return df, origin


def slice_image(image):
    indices = list(np.ndindex(*image.numblocks))
    slices = da.core.slices_from_chunks(image.chunks)
    img = [image[s] for s in slices]
    out = list(zip(img, indices))
    return out, slices


def slice_spots(spots, slices):
    bbox = [[[el.start, el.stop] for el in d] for d in slices]
    return [filter_spots(spots, b) for b in bbox]


def slice_data(image, spots):
    img_chunks, slices = slice_image(image)
    spot_chunks = slice_spots(spots, slices)
    return img_chunks, spot_chunks


def fit_old(spots, coo, scRNAseq, opts) -> Tuple[pd.DataFrame, pd.DataFrame]:
    label_image = coo.toarray()
    label_image = da.asarray(label_image)
    boundary = "reflect"
    label_image = da.overlap.overlap(label_image, 10, boundary)
    origin = spots[1]
    _spots = spots[0]
    out_1, out_2 = dask.delayed(fit_chunk, nout=2)(
        spots=_spots,
        coo=coo_matrix(label_image),
        scRNAseq=scRNAseq,
        opts=opts)
    cellData = da.from_delayed(out_1, shape=(30000, 12), dtype=object)
    # geneData = da.from_delayed(out_2)
    return cellData


def fit(spots, image, scRNAseq, opts):
    spots = dd.from_pandas(spots, npartitions=1)
    depth = 10
    image = da.asarray(image)
    # image = image[..., None]
    image = image.rechunk({1: 2000})
    image = image.rechunk({0: -1})

    boundary = "reflect"
    image = da.overlap.overlap(image, depth=depth, boundary=boundary)

    out_1,  out_2 = slice_data(image, spots)
    celltyped_blocks = np.empty(image.numblocks, dtype=object)

    lazy_data = []
    for img_block, spot_block in zip(out_1, out_2):
        block_id = img_block[1]
        spots_origin = spot_block[1]
        so = da.asarray(spots_origin)
        print(block_id)
        cd, gd = dask.delayed(fit_chunk, nout=2)(
            spots=spot_block[0],
            coo=coo_matrix(img_block[0]),
            scRNAseq=scRNAseq,
            opts=opts
        )
        # stick the spots origin next to the gene data inside a tuple.
        # Then get the cell data too
        lazy_data.append((cd, (gd, so)))

    return lazy_data






if __name__ == "__main__":
    from dask.distributed import Client, progress
    from pciSeq.src.core.logger import attach_to_log

    ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

    # set up the logger
    attach_to_log()

    # read some demo data
    _iss_spots = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'iss', 'spots.csv'))

    _coo = load_npz(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'segmentation', 'label_image.coo.npz'))

    _scRNAseq = pd.read_csv(os.path.join(ROOT_DIR, 'data', 'mouse', 'ca1', 'scRNA', 'scRNAseq.csv.gz'),
                            header=None, index_col=0, compression='gzip', dtype=object)
    _scRNAseq = _scRNAseq.rename(columns=_scRNAseq.iloc[0], copy=False).iloc[1:]
    _scRNAseq = _scRNAseq.astype(float).astype(np.uint32)

    # main task
    # _opts = {'max_iter': 10}
    _opts = {'save_data': True,
             'launch_viewer': False,
             'launch_diagnostics': False,
             # 'InsideCellBonus': 0,
             # 'cell_radius': 8.5
             }

    label_image = _coo.toarray()

    # client = Client(threads_per_worker=4, n_workers=1)
    tic = time.time()
    lazy_results = fit(_iss_spots, label_image, _scRNAseq, _opts)
    # client = Client(processes=False, n_workers=1, threads_per_worker=1)
    # print(client)
    with dask.config.set(scheduler='threads'):
        res = dask.compute(*lazy_results)
    toc = time.time()
    print(f"Computation time: {toc - tic:.2f} seconds\n")


    # label_image = da.overlap.overlap(da.asarray(label_image), 10, "reflect")
    # img, spt = chunk_data(label_image, _iss_spots)
    #
    # for i, input_block in enumerate(zip(spt, img)):
    #     _spots = input_block[0]
    #     _img = input_block[1]
    #     out = fit(_spots, _img, _scRNAseq, opts=_opts)
    #     print('done')
    # res = out.compute()


    print('ok')
    # fit(spots=_iss_spots, coo=_coo, scRNAseq=_scRNAseq, opts=_opts)

