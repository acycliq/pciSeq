"""Export a pciSeq run as a SpatialData object.

The store holds the spots as a points element, the segmentation as a labels
element, the cell typing results as an AnnData table annotating those labels,
and the scRNAseq reference as a second table. Run provenance (version, commit,
config) goes into the attrs dict.

Everything is written in its native acquisition space: x and y in pixels, z as
the plane index. One extra coordinate system called "microns" carries a Scale
transformation built from voxel_size, so the elements line up in physical
units without anybody having to know how pciSeq scales its coordinates
internally.

Two things are deliberately not in here. The mbtiles pyramid, because it is a
lossy display cache nothing in the scverse world can read, and it is bigger
than the raw image it came from. And the heavy per-cell diagnostics
(scaled_means and friends), which would roughly double the store to save a
multiplication; they stay in diagnostics.db.
"""
import json
import logging
import os
import shutil
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix

# spatialdata, anndata, dask and geopandas are imported inside the functions
# on purpose: core.io.export pulls this module in at import time, so hoisting them
# would make every `import pciSeq` pay for the whole scverse stack.
logger = logging.getLogger(__name__)


def _scale_to_microns(voxel_size, axes):
    """A transformation from native space (pixels and plane index) to microns.

    Native x and y are pixel positions, so they scale by the pixel size. Native
    z is the plane index, so it scales by the plane spacing. That is the whole
    anisotropy story as far as the store is concerned: the data stays as it was
    acquired and this transformation says what a step of one unit means.
    """
    from spatialdata.transformations import Scale

    per_axis = {'x': voxel_size[0], 'y': voxel_size[1], 'z': voxel_size[2]}
    return Scale([per_axis[ax] for ax in axes], axes=axes)


def _points(geneData: pd.DataFrame, varBayes, voxel_size):
    """The spots as a points element.

    Coordinates are x, y in pixels and the plane index as z. geneData also has
    a 'z' column but that one is the anisotropy scaled value the model works
    with, so it stays out of the store.

    The ragged per-spot lists (neighbour_array, neighbour_prob) are dropped:
    parquet handles them badly and the useful scalar, the probability of the
    assigned cell, is computed here instead.
    """
    from spatialdata.models import PointsModel

    is3D = 'z' in geneData.columns

    df = pd.DataFrame({
        'x': geneData.x.values.astype(np.float32),
        'y': geneData.y.values.astype(np.float32),
        'gene_name': geneData.gene_name.values,
        'spot_id': geneData.spot_id.values,
        'neighbour': geneData.neighbour.values,
        # the cell the spot physically sits in, 0 for none. Not the same as
        # 'neighbour', which is the cell the model picked.
        'inside_cell': geneData.inside_cell.values,
    })
    if is3D:
        df['z'] = geneData.plane_id.values.astype(np.float32)
        df['omp_score'] = geneData.omp_score.values.astype(np.float32)
        df['omp_intensity'] = geneData.omp_intensity.values.astype(np.float32)

    # probability of the assigned cell, and whether that cell is the
    # background. geneData rows are in the same order as the spots object.
    prob = varBayes.spots.parent_cell_prob
    df['neighbour_prob'] = prob.max(axis=1).astype(np.float32)
    df['is_hard_misread'] = (np.argmax(prob, axis=1) == prob.shape[1] - 1)

    axes = ('x', 'y', 'z') if is3D else ('x', 'y')
    coords = {ax: ax for ax in axes}
    return PointsModel.parse(
        df,
        coordinates=coords,
        feature_key='gene_name',
        transformations={'microns': _scale_to_microns(voxel_size, axes)},
    )


def _labels(coo: List[coo_matrix], label_map: Optional[Dict], voxel_size):
    """The segmentation masks as a labels element.

    Inside fit() the coo list holds the labels as pciSeq renumbered them, but the
    rest of the outputs were mapped back to the user's original labels before
    saving, so the masks are mapped back too. Otherwise the table would say cell
    2627 and the pixels would say 2176 and nothing would line up.

    Somebody calling this themselves may well hand over their own segmentation,
    with the original labels still in it. Mapping that back a second time would
    scramble it, so which of the two it is gets worked out first, see below.

    Built as a dask array one plane at a time because a full stack does not
    always fit in memory (a 6431 x 8544 x 68 uint32 stack is 15 GB dense).
    Zarr then writes it plane by plane as well.
    """
    import dask
    import dask.array as da
    import fastremap
    from spatialdata.models import Labels2DModel, Labels3DModel

    # label_map goes original -> renumbered, we need the other direction. But only if
    # the masks really are in the renumbered labels. pciSeq renumbers only when the
    # labels are not sequential, and then the n cells become exactly 1..n. n labels
    # that are not 1..n must have one above n, so the biggest label settles it: up to
    # n means renumbered, above n means these are the originals already.
    inverse = None
    if label_map:
        n_cells = sum(1 for original in label_map if original != 0)
        biggest = max((int(m.data.max()) for m in coo if m.nnz), default=0)
        if biggest <= n_cells:
            inverse = {v: k for k, v in label_map.items()}
        else:
            logger.info('the masks already hold the original labels, not mapping them back')

    def plane(sparse_plane):
        arr = sparse_plane.toarray().astype(np.uint32)
        if inverse:
            fastremap.remap(arr, inverse, preserve_missing_labels=True, in_place=True)
        return arr

    shape = coo[0].shape
    planes = [
        da.from_delayed(dask.delayed(plane)(p), shape=shape, dtype=np.uint32)
        for p in coo
    ]

    if len(planes) == 1:
        return Labels2DModel.parse(
            planes[0], dims=('y', 'x'),
            transformations={'microns': _scale_to_microns(voxel_size, ('y', 'x'))},
        )
    return Labels3DModel.parse(
        da.stack(planes, axis=0), dims=('z', 'y', 'x'),
        transformations={'microns': _scale_to_microns(voxel_size, ('z', 'y', 'x'))},
    )


def _counts_matrix(cellData: pd.DataFrame, gene_panel) -> csr_matrix:
    """The cell by gene counts, unpacked from the ragged Genenames and
    CellGeneCount lists into a proper sparse matrix with one column per gene
    in the panel."""
    col_of = {g: j for j, g in enumerate(gene_panel)}
    rows, cols, vals = [], [], []
    for i, (names, counts) in enumerate(zip(cellData.Genenames, cellData.CellGeneCount)):
        for g, c in zip(names, counts):
            rows.append(i)
            cols.append(col_of[g])
            vals.append(c)
    return coo_matrix(
        (np.asarray(vals, dtype=np.float32), (rows, cols)),
        shape=(len(cellData), len(gene_panel)),
    ).tocsr()


def _cells_table(cellData: pd.DataFrame, varBayes, cfg, labels_name: str):
    """The cell typing results as an AnnData table annotating the labels.

    X is the sparse cell by gene count matrix, var is the gene panel, obs the
    headline results, obsm the per-cell arrays. classProb row 0 is the
    background pseudocell, which cellData does not carry, so the model arrays
    are sliced from 1 to line up with the cellData rows.
    """
    from anndata import AnnData
    from spatialdata.models import TableModel

    n = len(cellData)
    gene_panel = varBayes.genes.gene_panel
    class_names = np.asarray(varBayes.cells.class_names, dtype=str)
    classProb = varBayes.cells.classProb[1:n + 1]

    top = classProb.argmax(axis=1)
    obs = pd.DataFrame({
        'cell_num': cellData.Cell_Num.values,
        'region': pd.Categorical([labels_name] * n),
        'class_name': class_names[top],
        'class_prob': classProb[np.arange(n), top].astype(np.float32),
    })
    obs.index = cellData.Cell_Num.astype(str).values

    # the tie freezing state: which cells had their class pinned, and when
    freezer = getattr(varBayes, 'tie_freezer', None)
    if freezer is not None:
        obs['is_pinned'] = freezer.frozen[1:n + 1]
        pinned_at = np.full(n, -1, dtype=np.int32)
        for c, snap in freezer.snapshots.items():
            if 1 <= c <= n:
                pinned_at[c - 1] = snap['iteration']
        obs['pinned_at_iteration'] = pinned_at

    # centroids in native space. cellData X, Y, Z are anisotropy scaled, which
    # only differs from pixels on the z axis while the pixels are square, so z
    # is divided back into plane units.
    vs = cfg['voxel_size']
    if 'Z' in cellData.columns:
        z_native = cellData.Z.values / (vs[2] / vs[0])
        spatial = np.column_stack([cellData.X.values, cellData.Y.values, z_native])
    else:
        spatial = np.column_stack([cellData.X.values, cellData.Y.values])

    adata = AnnData(
        X=_counts_matrix(cellData, gene_panel),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(gene_panel, name='gene_name')),
    )
    adata.obsm['spatial'] = spatial.astype(np.float32)
    adata.obsm['class_prob'] = classProb.astype(np.float32)
    adata.uns['class_names'] = class_names.tolist()

    return TableModel.parse(
        adata, region=labels_name, region_key='region', instance_key='cell_num',
    )


def _reference_table(varBayes):
    """The scRNAseq reference the run was called against, classes as rows and
    genes as columns, so the store says what produced these calls."""
    from anndata import AnnData

    ref = varBayes.single_cell.raw_data  # genes x classes
    return AnnData(
        X=np.asarray(ref.values, dtype=np.float32).T,
        obs=pd.DataFrame(index=pd.Index(ref.columns.astype(str), name='class_name')),
        var=pd.DataFrame(index=pd.Index(ref.index.astype(str), name='gene_name')),
    )


def _attrs(varBayes, cfg) -> Dict:
    """Provenance: enough to trace the store back to the code and config that
    produced it. label_map is dropped because it can hold one entry per cell
    and the store already carries the original labels everywhere."""
    import spatialdata

    config = {}
    for k, v in cfg.items():
        if k == 'label_map':
            continue
        try:
            json.dumps(v)
            config[k] = v
        except TypeError:
            config[k] = str(v)

    return {
        'pciseq': {
            **varBayes.metadata,
            'config': config,
            'spatialdata_version': spatialdata.__version__,
        }
    }


def to_spatialdata(cellData: pd.DataFrame,
                   geneData: pd.DataFrame,
                   coo: List[coo_matrix],
                   varBayes,
                   cfg: Dict):
    """Assemble the run into an in-memory SpatialData object.

    Args:
        cellData: cell typing results, one row per cell, original labels.
        geneData: spot results, one row per spot, original labels.
        coo: the segmentation, one sparse plane per z. Either the list fit()
            renumbered in place or a copy with the original labels; which of the
            two it is gets detected, and the store always ends up with the
            original labels.
        varBayes: the fitted model, read for the arrays the two frames do not
            carry (class posterior, spot probabilities, gene panel, reference).
        cfg: the resolved config. voxel_size and label_map are used here.
    """
    from spatialdata import SpatialData

    voxel_size = cfg['voxel_size']
    labels_name = 'cell_labels'

    sdata = SpatialData(
        points={'transcripts': _points(geneData, varBayes, voxel_size)},
        labels={labels_name: _labels(coo, cfg.get('label_map'), voxel_size)},
        tables={
            'cells': _cells_table(cellData, varBayes, cfg, labels_name),
            'reference': _reference_table(varBayes),
        },
    )
    sdata.attrs = _attrs(varBayes, cfg)
    return sdata


def _store_voxel_size(sdata, voxel_size):
    """The store's own voxel size, or the caller's override.

    Reading it from the provenance keeps a late insert consistent with the
    elements already in there: the microns transform cannot contradict the one
    the run was written with.
    """
    if voxel_size is not None:
        return voxel_size
    try:
        return sdata.attrs['pciseq']['config']['voxel_size']
    except KeyError:
        raise ValueError('the store carries no voxel_size, pass voxel_size=[x, y, z]')


def add_image(store_path: str,
              img,
              name: str = 'background',
              voxel_size=None,
              scale_factors=(2, 2, 2)) -> None:
    """Add a background image to an existing SpatialData store.

    Works like inserting into a database: the store already holds the run
    (spots, labels, tables) and this writes one more element into it, without
    touching anything else. Kept separate from write_spatialdata because fit()
    never sees the image, the same reason stage_image is its own entry point.

    Args:
        store_path: path to an existing spatialdata.zarr written by pciSeq.
        img: the image to add. A 2D array (H, W), a 3D stack (Z, H, W), or a
            3D stack with channels (Z, H, W, C). Same shapes stage_image takes.
        name: element name inside the store.
        voxel_size: [x, y, z]. Left as None it is read from the store's own
            provenance, so the microns transform automatically matches the
            other elements. Pass it only for a store that lacks the metadata.
        scale_factors: the multiscale pyramid, each level relative to the one
            before. The default (2, 2, 2) gives four scales. None writes a
            single scale.
    """
    from spatialdata import read_zarr
    from spatialdata.models import Image2DModel, Image3DModel

    sdata = read_zarr(store_path)
    voxel_size = _store_voxel_size(sdata, voxel_size)

    img = np.asarray(img)
    if img.ndim == 2:
        # (H, W) -> (c, y, x)
        parsed = Image2DModel.parse(
            img[None, :, :], dims=('c', 'y', 'x'),
            scale_factors=list(scale_factors) if scale_factors else None,
            transformations={'microns': _scale_to_microns(voxel_size, ('y', 'x'))},
        )
    else:
        if img.ndim == 3:
            data = img[None, :, :, :]              # (Z, H, W)   -> (c, z, y, x)
        elif img.ndim == 4:
            data = np.moveaxis(img, -1, 0)         # (Z, H, W, C) -> (c, z, y, x)
        else:
            raise ValueError(f'expected 2, 3 or 4 dimensions, got {img.ndim}')
        # downsample y and x only. z is always the thin axis (tens of planes
        # against thousands of pixels) and halving it as well runs out of
        # planes after a couple of levels.
        factors = ([{'z': 1, 'y': f, 'x': f} for f in scale_factors]
                   if scale_factors else None)
        parsed = Image3DModel.parse(
            data, dims=('c', 'z', 'y', 'x'),
            scale_factors=factors,
            transformations={'microns': _scale_to_microns(voxel_size, ('z', 'y', 'x'))},
        )

    sdata.images[name] = parsed
    sdata.write_element(name)
    logger.info('image %r added to %s', name, store_path)


def add_boundaries(store_path: str,
                   labels_name: str = 'cell_labels',
                   name: str = 'cell_boundaries',
                   voxel_size=None) -> None:
    """Add per-plane cell boundary polygons to an existing SpatialData store.

    Same insert pattern as add_image, but this one needs no data at all: the
    boundaries are derived from the segmentation, and the store already holds
    the segmentation as the labels element. They are extracted here with the
    same chain code tracing the pipeline uses, so what goes in matches what
    pciSeq would have drawn.

    Shapes in SpatialData are strictly 2D, so a 3D run gets one shapes element
    per plane, cell_boundaries_plane_000 and so on, each a set of polygons
    indexed by the cell label. A cell spanning 12 planes appears as 12
    polygons, one per element, all under its own label, which matches the
    per-plane boundary files pciSeq writes. A 2D run gets a single element.

    Args:
        store_path: path to an existing spatialdata.zarr written by pciSeq.
        labels_name: the labels element to trace.
        name: element name, used as a prefix on 3D runs.
        voxel_size: [x, y, z]. Left as None it is read from the store's own
            provenance, same as add_image.
    """
    import geopandas as gpd
    from shapely.geometry import Polygon
    from spatialdata import read_zarr
    from spatialdata.models import ShapesModel

    from ...preprocess.cell_processing import extract_borders

    sdata = read_zarr(store_path)

    labels = sdata.labels[labels_name]
    if hasattr(labels, 'scale0'):  # multiscale, take the full resolution
        labels = next(iter(labels['scale0'].values()))
    stack = np.asarray(labels)
    if stack.ndim == 2:
        stack = stack[None, :, :]

    voxel_size = _store_voxel_size(sdata, voxel_size)
    is3D = stack.shape[0] > 1

    written = []
    for z in range(stack.shape[0]):
        borders = extract_borders(stack[z].astype(np.uint32))
        polys, labs = [], []
        for lab, coords in zip(borders.label, borders.coords):
            if len(coords) < 4:
                # a closed ring needs at least three distinct points, cells of
                # one or two pixels on this plane have no polygon to draw
                continue
            polys.append(Polygon(coords))
            labs.append(int(lab))
        if not polys:
            continue

        gdf = gpd.GeoDataFrame({'geometry': polys}, index=pd.Index(labs, name='label'))
        element_name = f'{name}_plane_{z:03d}' if is3D else name
        sdata.shapes[element_name] = ShapesModel.parse(
            gdf,
            transformations={'microns': _scale_to_microns(voxel_size, ('x', 'y'))},
        )
        sdata.write_element(element_name)
        written.append(element_name)

    logger.info('%d boundary element(s) added to %s', len(written), store_path)


def write_spatialdata(cellData, geneData, coo, varBayes, cfg, out_dir: str) -> str:
    """Write the run as a SpatialData zarr store and return its path."""
    path = os.path.join(out_dir, 'spatialdata.zarr')
    if os.path.exists(path):
        shutil.rmtree(path)

    sdata = to_spatialdata(cellData, geneData, coo, varBayes, cfg)
    sdata.write(path)
    logger.info('SpatialData store saved at %s', path)
    return path
