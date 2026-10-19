"""Tests for the SpatialData export.

Built around one small synthetic run: three cells on a 64 x 64 grid over two
planes, a handful of spots, a three class reference. The point is not the
biology, it is that every piece of the store lines up with every other piece:
the table's cell numbers match the pixel values in the labels, the spots land
inside the cells they were assigned to, and the microns coordinate system
scales everything consistently.
"""

import types

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import coo_matrix

from pciSeq.src.core.utils.spatialdata_export import (
    add_boundaries, add_image, to_spatialdata, write_spatialdata)

VOXEL = [0.5, 0.5, 2.0]  # x, y, z. Sz = 4, deliberately not 1


def _synthetic_run():
    """A tiny fake run, shaped exactly like what fit() has at the save step.

    The segmentation labels are 10, 20, 30 (deliberately non sequential) and
    label_map records how pciSeq renumbered them to 1, 2, 3. The coo planes
    hold the renumbered labels, cellData and geneData the original ones, which
    is the state the exporter really receives.
    """
    label_map = {10: 1, 20: 2, 30: 3}

    # three 8x8 square cells along the diagonal, same on both planes
    plane = np.zeros((64, 64), dtype=np.uint32)
    for lab, (r, c) in zip([1, 2, 3], [(8, 8), (28, 28), (48, 48)]):
        plane[r:r + 8, c:c + 8] = lab
    coo = [coo_matrix(plane), coo_matrix(plane)]

    # spots sitting in the middle of their cells (x is the column axis)
    geneData = pd.DataFrame({
        'gene_name': ['Vip', 'Sst', 'Vip', 'Pvalb'],
        'gene_id': [2, 1, 2, 0],
        'spot_id': [0, 1, 2, 3],
        'x': np.array([12.0, 32.0, 33.0, 52.0]),
        'y': np.array([12.0, 32.0, 31.0, 52.0]),
        'z': np.array([0.0, 4.0, 4.0, 4.0]),  # the scaled z the model uses
        'plane_id': np.array([0, 1, 1, 1]),
        'neighbour': [10, 20, 20, 30],
        # the cell each spot physically sits in. Same as neighbour here because the
        # spots were put in the middle of their cells, but it is a different thing.
        'inside_cell': [10, 20, 20, 30],
        'omp_score': [0.9, 0.8, 0.7, 0.6],
        'omp_intensity': [1.0, 2.0, 3.0, 4.0],
    })

    # cellData as collect_data makes it: ragged lists, original labels,
    # centroids in the scaled space (z multiplied by Sz = 4)
    cellData = pd.DataFrame({
        'Cell_Num': [10, 20, 30],
        'X': [12.0, 32.0, 52.0],
        'Y': [12.0, 32.0, 52.0],
        'Z': [2.0, 2.0, 2.0],  # plane 0.5, times Sz
        'Genenames': [['Vip'], ['Sst', 'Vip'], ['Pvalb']],
        'CellGeneCount': [[1.0], [1.0, 1.0], [1.0]],
        'ClassName': [['A'], ['B'], ['Zero']],
        'Prob': [[0.9], [0.8], [0.7]],
    })

    # the model object, only the bits the exporter reads. classProb row 0 is
    # the background pseudocell.
    classProb = np.array([
        [0.0, 0.0, 1.0],
        [0.9, 0.05, 0.05],
        [0.1, 0.8, 0.1],
        [0.2, 0.1, 0.7],
    ], dtype=np.float32)
    parent_cell_prob = np.array([
        [0.9, 0.1],
        [0.8, 0.2],
        [0.7, 0.3],
        [0.1, 0.9],  # this spot's best parent is the background column
    ], dtype=np.float32)

    ref = pd.DataFrame(
        np.arange(9, dtype=np.float32).reshape(3, 3),
        index=pd.Index(['Pvalb', 'Sst', 'Vip'], name='gene_name'),
        columns=['A', 'B', 'Zero'],
    )

    varBayes = types.SimpleNamespace(
        cells=types.SimpleNamespace(
            classProb=classProb,
            class_names=np.array(['A', 'B', 'Zero']),
        ),
        spots=types.SimpleNamespace(parent_cell_prob=parent_cell_prob),
        genes=types.SimpleNamespace(gene_panel=np.array(['Pvalb', 'Sst', 'Vip'])),
        single_cell=types.SimpleNamespace(raw_data=ref),
        metadata={'version': '0.0.66.dev3', 'branch': 'test', 'commit': 'abc123',
                  'commit_date': 'today', 'created_at': 'now'},
    )

    cfg = {'voxel_size': VOXEL, 'label_map': label_map, 'is3D': True,
           'save_data': True, 'output_path': 'default'}
    return cellData, geneData, coo, varBayes, cfg


@pytest.fixture(scope='module')
def sdata():
    return to_spatialdata(*_synthetic_run())


class TestElements:

    def test_the_store_has_all_four_elements(self, sdata):
        assert 'transcripts' in sdata.points
        assert 'cell_labels' in sdata.labels
        assert 'cells' in sdata.tables
        assert 'reference' in sdata.tables

    def test_labels_carry_the_original_segmentation_labels(self, sdata):
        # the coo planes hold the renumbered 1, 2, 3 and the exporter must map
        # them back, otherwise nothing matches the table
        vals = set(np.unique(np.asarray(sdata.labels['cell_labels'])))
        assert vals == {0, 10, 20, 30}

    def test_masks_with_the_original_labels_are_not_mapped_back_again(self):
        # somebody calling the exporter themselves hands over their own segmentation
        # with the original labels still in it. The originals here are 2, 3, 7 and
        # pciSeq renumbered them to 1, 2, 3, so an original 2 is also a renumbered 2.
        # Mapping back a second time read it as the renumbered one and made it a 3,
        # and made the 3 a 7: two cells relabelled and the table no longer matching.
        from pciSeq.src.core.utils.spatialdata_export import _labels
        plane = np.zeros((16, 16), dtype=np.uint32)
        plane[2:5, 2:5], plane[7:10, 7:10], plane[12:15, 12:15] = 2, 3, 7
        label_map = {0: 0, 2: 1, 3: 2, 7: 3}

        out = np.asarray(_labels([coo_matrix(plane), coo_matrix(plane)], label_map, VOXEL))
        assert np.array_equal(out[0], plane)
        assert np.array_equal(out[1], plane)

    def test_the_renumbered_masks_of_that_same_run_are_mapped_back(self):
        # the other way round, which is what fit() hands over: same cells, labels 1, 2, 3
        from pciSeq.src.core.utils.spatialdata_export import _labels
        original = np.zeros((16, 16), dtype=np.uint32)
        original[2:5, 2:5], original[7:10, 7:10], original[12:15, 12:15] = 2, 3, 7
        renumbered = np.zeros((16, 16), dtype=np.uint32)
        renumbered[2:5, 2:5], renumbered[7:10, 7:10], renumbered[12:15, 12:15] = 1, 2, 3
        label_map = {0: 0, 2: 1, 3: 2, 7: 3}

        out = np.asarray(_labels([coo_matrix(renumbered)], label_map, VOXEL))
        assert np.array_equal(out, original)

    def test_without_a_label_map_nothing_is_touched(self):
        from pciSeq.src.core.utils.spatialdata_export import _labels
        plane = np.zeros((16, 16), dtype=np.uint32)
        plane[2:5, 2:5], plane[7:10, 7:10] = 1, 2
        out = np.asarray(_labels([coo_matrix(plane)], None, VOXEL))
        assert np.array_equal(out, plane)

    def test_points_use_the_plane_index_not_the_scaled_z(self, sdata):
        df = sdata.points['transcripts'].compute()
        assert sorted(df.z.unique()) == [0.0, 1.0]

    def test_hard_misread_flag(self, sdata):
        df = sdata.points['transcripts'].compute()
        assert df.is_hard_misread.tolist() == [False, False, False, True]

    def test_counts_matrix_is_cells_by_panel(self, sdata):
        adata = sdata.tables['cells']
        assert adata.shape == (3, 3)
        assert adata.var_names.tolist() == ['Pvalb', 'Sst', 'Vip']
        # cell 20 (row 1) had Sst and Vip, one read each
        assert adata.X[1].toarray().tolist() == [[0.0, 1.0, 1.0]]

    def test_table_annotates_the_labels(self, sdata):
        adata = sdata.tables['cells']
        assert adata.obs['cell_num'].tolist() == [10, 20, 30]
        attrs = adata.uns['spatialdata_attrs']
        assert attrs['region'] == 'cell_labels'
        assert attrs['instance_key'] == 'cell_num'

    def test_class_posterior_and_headline_agree(self, sdata):
        adata = sdata.tables['cells']
        assert adata.obs['class_name'].tolist() == ['A', 'B', 'Zero']
        assert adata.obsm['class_prob'].shape == (3, 3)
        got = adata.obsm['class_prob'].argmax(axis=1)
        names = np.array(adata.uns['class_names'])[got]
        assert names.tolist() == adata.obs['class_name'].tolist()

    def test_centroid_z_is_back_in_plane_units(self, sdata):
        spatial = sdata.tables['cells'].obsm['spatial']
        # cellData carried Z = 2.0 in scaled units, Sz = 4, so 0.5 planes
        assert spatial[:, 2] == pytest.approx([0.5, 0.5, 0.5])

    def test_provenance(self, sdata):
        meta = sdata.attrs['pciseq']
        assert meta['commit'] == 'abc123'
        assert 'label_map' not in meta['config']
        assert meta['config']['voxel_size'] == VOXEL


class TestConsistency:
    """The checks that would catch a wrong transform, which is the failure
    mode that looks plausible and is silently wrong."""

    def test_every_spot_lands_inside_its_assigned_cell(self, sdata):
        labels = np.asarray(sdata.labels['cell_labels'])  # (z, y, x)
        df = sdata.points['transcripts'].compute()
        for _, s in df[~df.is_hard_misread].iterrows():
            pixel = labels[int(round(s.z)), int(round(s.y)), int(round(s.x))]
            assert pixel == s.neighbour, \
                f"spot {int(s.spot_id)} assigned to {s.neighbour} but sits on pixel {pixel}"

    def test_microns_transform_scales_all_elements_consistently(self, sdata):
        from spatialdata.transformations import get_transformation

        for name, element in [('points', sdata.points['transcripts']),
                              ('labels', sdata.labels['cell_labels'])]:
            t = get_transformation(element, 'microns').to_affine_matrix(
                input_axes=('x', 'y', 'z'), output_axes=('x', 'y', 'z'))
            assert np.diag(t)[:3] == pytest.approx(VOXEL), name

    def test_spot_and_its_cell_centroid_agree_in_microns(self, sdata):
        # spot 0 sits at the centre of cell 10. Carrying both into microns
        # through their own routes must land them at the same place.
        df = sdata.points['transcripts'].compute()
        spot = df[df.spot_id == 0].iloc[0]
        spot_um = np.array([spot.x, spot.y, spot.z]) * np.array(VOXEL)

        adata = sdata.tables['cells']
        centroid = adata.obsm['spatial'][adata.obs.cell_num.tolist().index(10)]
        centroid_um = centroid * np.array(VOXEL)

        assert np.abs(spot_um - centroid_um).max() < max(VOXEL) * 2


class TestAddImage:
    """Adding a background image to an already written store, the database
    insert pattern: the run is saved first, the image arrives later."""

    @pytest.fixture()
    def store(self, tmp_path):
        cellData, geneData, coo, varBayes, cfg = _synthetic_run()
        return write_spatialdata(cellData, geneData, coo, varBayes, cfg, str(tmp_path))

    def test_image_lands_next_to_the_run(self, store):
        import spatialdata

        img = np.random.default_rng(0).integers(0, 255, (2, 64, 64), dtype=np.uint8)
        add_image(store, img, scale_factors=(2,))

        back = spatialdata.read_zarr(store)
        assert 'background' in back.images
        # nothing else was disturbed
        assert back.tables['cells'].shape == (3, 3)
        assert len(back.points['transcripts'].compute()) == 4

    def test_voxel_size_comes_from_the_store_itself(self, store):
        import spatialdata
        from spatialdata.transformations import get_transformation

        img = np.zeros((2, 64, 64), dtype=np.uint8)
        add_image(store, img, name='dapi', scale_factors=None)

        back = spatialdata.read_zarr(store)
        t = get_transformation(back.images['dapi'], 'microns').to_affine_matrix(
            input_axes=('x', 'y', 'z'), output_axes=('x', 'y', 'z'))
        # the transform must match the store's own voxel_size, nobody passed it
        assert np.diag(t)[:3] == pytest.approx(VOXEL)

    def test_pyramid_is_built(self, store):
        import spatialdata

        img = np.zeros((2, 64, 64), dtype=np.uint8)
        add_image(store, img, name='pyr', scale_factors=(2, 2))

        back = spatialdata.read_zarr(store)
        levels = [k for k in back.images['pyr'].keys() if k.startswith('scale')]
        assert len(levels) == 3, "two scale_factors should give three levels"

    def test_2d_image(self, store):
        import spatialdata

        add_image(store, np.zeros((64, 64), dtype=np.uint8), name='flat',
                  scale_factors=None)
        back = spatialdata.read_zarr(store)
        assert back.images['flat'].dims == ('c', 'y', 'x')


class TestAddBoundaries:
    """The boundaries need no input data: they are traced from the labels the
    store already holds, one shapes element per plane."""

    @pytest.fixture()
    def store(self, tmp_path):
        cellData, geneData, coo, varBayes, cfg = _synthetic_run()
        return write_spatialdata(cellData, geneData, coo, varBayes, cfg, str(tmp_path))

    def test_one_element_per_plane_indexed_by_label(self, store):
        import spatialdata

        add_boundaries(store)
        back = spatialdata.read_zarr(store)

        assert 'cell_boundaries_plane_000' in back.shapes
        assert 'cell_boundaries_plane_001' in back.shapes
        for z in (0, 1):
            gdf = back.shapes[f'cell_boundaries_plane_{z:03d}']
            assert sorted(gdf.index.tolist()) == [10, 20, 30], \
                "polygons must be indexed by the original segmentation labels"

    def test_each_polygon_contains_its_cell_centroid(self, store):
        import spatialdata
        from shapely.geometry import Point

        add_boundaries(store)
        back = spatialdata.read_zarr(store)
        gdf = back.shapes['cell_boundaries_plane_000']
        centres = {10: (12, 12), 20: (32, 32), 30: (52, 52)}
        for lab, (cx, cy) in centres.items():
            assert gdf.loc[lab].geometry.contains(Point(cx, cy)), \
                f"boundary of cell {lab} should contain its centre"

    def test_rest_of_the_store_untouched(self, store):
        import spatialdata

        add_boundaries(store)
        back = spatialdata.read_zarr(store)
        assert back.tables['cells'].shape == (3, 3)
        assert len(back.points['transcripts'].compute()) == 4


class TestRoundTrip:

    def test_write_and_read_back(self, tmp_path):
        cellData, geneData, coo, varBayes, cfg = _synthetic_run()
        path = write_spatialdata(cellData, geneData, coo, varBayes, cfg, str(tmp_path))

        import spatialdata
        back = spatialdata.read_zarr(path)
        assert back.tables['cells'].shape == (3, 3)
        assert back.tables['cells'].obs['cell_num'].tolist() == [10, 20, 30]
        assert set(np.unique(np.asarray(back.labels['cell_labels']))) == {0, 10, 20, 30}
        assert len(back.points['transcripts'].compute()) == 4
        assert back.attrs['pciseq']['commit'] == 'abc123'
