"""
Tests for process_labels (sparse fastremap approach) and
calculate_cell_properties (numba approach).
Verifies both produce identical results to the original dense approaches.
"""

import numpy as np
import pytest
from scipy.sparse import coo_matrix
import fastremap
import skimage.measure as skmeas
from pciSeq.src.preprocess.label_processing import process_labels
from pciSeq.src.preprocess.cell_processing import calculate_cell_properties


def _process_labels_dense(coo_list):
    """Original dense approach for comparison."""
    arr_3d = np.stack([d.toarray() for d in coo_list])
    unique_labels = np.unique(arr_3d[arr_3d > 0])
    is_sequential = len(unique_labels) == unique_labels.max()

    if not is_sequential:
        normalized, label_map = fastremap.renumber(arr_3d, in_place=False, preserve_zero=True)
    else:
        normalized, label_map = arr_3d, None

    return [coo_matrix(d) for d in normalized], label_map


def _make_coo_list(masks):
    """Convert a 3D numpy array to a list of coo_matrix."""
    return [coo_matrix(d) for d in masks]


class TestProcessLabels:
    """Test that sparse process_labels matches the dense approach."""

    def test_sequential_labels_no_renumbering(self):
        """When labels are already 1..N, no renumbering should happen."""
        masks = np.zeros((3, 10, 10), dtype=np.uint16)
        masks[0, 1:4, 1:4] = 1
        masks[1, 5:8, 5:8] = 2
        masks[2, 2:5, 2:5] = 3

        coo_dense = _make_coo_list(masks.copy())
        coo_sparse = _make_coo_list(masks.copy())

        result_dense, map_dense = _process_labels_dense(coo_dense)
        result_sparse, map_sparse = process_labels(coo_sparse)

        assert map_dense is None
        assert map_sparse is None
        for d, s in zip(result_dense, result_sparse):
            np.testing.assert_array_equal(d.toarray(), s.toarray())

    def test_gapped_labels_renumbered(self):
        """Labels with gaps (e.g. 1, 5, 10) should be renumbered to 1..N."""
        masks = np.zeros((3, 10, 10), dtype=np.uint16)
        masks[0, 1:4, 1:4] = 5
        masks[1, 5:8, 5:8] = 10
        masks[2, 2:5, 2:5] = 20

        coo_dense = _make_coo_list(masks.copy())
        coo_sparse = _make_coo_list(masks.copy())

        result_dense, map_dense = _process_labels_dense(coo_dense)
        result_sparse, map_sparse = process_labels(coo_sparse)

        assert map_dense is not None
        assert map_sparse is not None
        for d, s in zip(result_dense, result_sparse):
            np.testing.assert_array_equal(d.toarray(), s.toarray())

    def test_shared_labels_across_planes(self):
        """Same label appearing in multiple planes (typical for 3D cells)."""
        masks = np.zeros((5, 20, 20), dtype=np.uint16)
        # Cell 3 spans planes 0-2
        masks[0, 2:6, 2:6] = 3
        masks[1, 2:6, 2:6] = 3
        masks[2, 3:5, 3:5] = 3
        # Cell 7 spans planes 2-4
        masks[2, 10:15, 10:15] = 7
        masks[3, 10:15, 10:15] = 7
        masks[4, 11:14, 11:14] = 7

        coo_dense = _make_coo_list(masks.copy())
        coo_sparse = _make_coo_list(masks.copy())

        result_dense, _ = _process_labels_dense(coo_dense)
        result_sparse, _ = process_labels(coo_sparse)

        for d, s in zip(result_dense, result_sparse):
            np.testing.assert_array_equal(d.toarray(), s.toarray())

    def test_empty_planes(self):
        """Planes with no cells should be handled gracefully."""
        masks = np.zeros((5, 10, 10), dtype=np.uint16)
        # Only planes 1 and 3 have cells
        masks[1, 2:5, 2:5] = 4
        masks[3, 6:9, 6:9] = 8

        coo_dense = _make_coo_list(masks.copy())
        coo_sparse = _make_coo_list(masks.copy())

        result_dense, _ = _process_labels_dense(coo_dense)
        result_sparse, _ = process_labels(coo_sparse)

        for d, s in zip(result_dense, result_sparse):
            np.testing.assert_array_equal(d.toarray(), s.toarray())

    def test_many_labels_with_gaps(self):
        """Realistic scenario: many labels with random gaps."""
        rng = np.random.default_rng(42)
        masks = np.zeros((10, 50, 50), dtype=np.uint16)

        # Place 50 cells with random non-sequential labels
        labels = sorted(rng.choice(range(1, 500), size=50, replace=False))
        for i, label in enumerate(labels):
            plane = i % 10
            r, c = rng.integers(0, 40, size=2)
            masks[plane, r:r+5, c:c+5] = label

        coo_dense = _make_coo_list(masks.copy())
        coo_sparse = _make_coo_list(masks.copy())

        result_dense, map_dense = _process_labels_dense(coo_dense)
        result_sparse, map_sparse = process_labels(coo_sparse)

        assert map_dense is not None
        assert map_sparse is not None
        for d, s in zip(result_dense, result_sparse):
            np.testing.assert_array_equal(d.toarray(), s.toarray())

    def test_single_plane(self):
        """Edge case: only one plane."""
        masks = np.zeros((1, 10, 10), dtype=np.uint16)
        masks[0, 1:4, 1:4] = 5
        masks[0, 6:9, 6:9] = 12

        coo_dense = _make_coo_list(masks.copy())
        coo_sparse = _make_coo_list(masks.copy())

        result_dense, _ = _process_labels_dense(coo_dense)
        result_sparse, _ = process_labels(coo_sparse)

        for d, s in zip(result_dense, result_sparse):
            np.testing.assert_array_equal(d.toarray(), s.toarray())

    def test_background_preserved(self):
        """Background (label 0) should never be remapped."""
        masks = np.zeros((3, 10, 10), dtype=np.uint16)
        masks[0, 1:3, 1:3] = 10
        masks[1, 5:7, 5:7] = 20

        _, _ = process_labels(_make_coo_list(masks.copy()))

        # After renumbering, the sparse matrices should have no zeros in .data
        coo_list = _make_coo_list(masks.copy())
        result, _ = process_labels(coo_list)
        for m in result:
            if m.nnz > 0:
                assert np.all(m.data > 0), "Background label 0 should not appear in sparse .data"


def _calculate_cell_properties_regionprops(masks, voxel_size):
    """Original regionprops approach for comparison."""
    scaling = [voxel_size[0] / voxel_size[0], voxel_size[1] / voxel_size[0], voxel_size[2] / voxel_size[0]]
    scaling = scaling[::-1]

    props = skmeas.regionprops_table(
        label_image=masks,
        spacing=scaling,
        properties=['label', 'area', 'centroid', 'bbox']
    )

    import pandas as pd
    props_df = pd.DataFrame(props)
    props_df['mean_area_per_slice'] = (
        props_df['area'].values /
        (props_df['bbox-3'].values - props_df['bbox-0'].values)
    )
    props_df = props_df.rename(columns={
        "mean_area_per_slice": 'area', 'area': 'volume',
        'centroid-0': 'z_cell', 'centroid-1': 'y_cell', 'centroid-2': 'x_cell'
    })
    props_df = props_df[['label', 'area', 'z_cell', 'y_cell', 'x_cell']]
    return props_df.astype({
        "label": np.uint32, "area": np.uint32,
        'z_cell': np.float32, 'y_cell': np.float32, 'x_cell': np.float32
    })


class TestCellProperties:
    """Test that numba-based calculate_cell_properties matches regionprops."""

    def test_basic_cells(self):
        """Simple cells across multiple planes."""
        masks = np.zeros((5, 30, 30), dtype=np.uint32)
        masks[0, 2:8, 2:8] = 1
        masks[1, 2:8, 2:8] = 1
        masks[2, 3:7, 3:7] = 1
        masks[2, 15:25, 15:25] = 2
        masks[3, 15:25, 15:25] = 2
        masks[4, 16:24, 16:24] = 2

        voxel_size = [0.28, 0.28, 0.7]
        expected = _calculate_cell_properties_regionprops(masks, voxel_size)
        result = calculate_cell_properties(_make_coo_list(masks), voxel_size)

        assert len(result) == len(expected)
        for col in ['label', 'area', 'x_cell', 'y_cell', 'z_cell']:
            np.testing.assert_array_equal(result[col].values, expected[col].values)

    def test_isotropic_voxels(self):
        """Voxel size [1, 1, 1] — no scaling."""
        masks = np.zeros((3, 20, 20), dtype=np.uint32)
        masks[0, 5:10, 5:10] = 1
        masks[1, 5:10, 5:10] = 1
        masks[2, 5:10, 5:10] = 1

        voxel_size = [1, 1, 1]
        expected = _calculate_cell_properties_regionprops(masks, voxel_size)
        result = calculate_cell_properties(_make_coo_list(masks), voxel_size)

        for col in ['label', 'area', 'x_cell', 'y_cell', 'z_cell']:
            np.testing.assert_array_equal(result[col].values, expected[col].values)

    def test_single_plane_cell(self):
        """Cell existing in only one plane."""
        masks = np.zeros((5, 20, 20), dtype=np.uint32)
        masks[2, 3:8, 3:8] = 1
        masks[0, 12:18, 12:18] = 2
        masks[1, 12:18, 12:18] = 2

        voxel_size = [0.28, 0.28, 0.7]
        expected = _calculate_cell_properties_regionprops(masks, voxel_size)
        result = calculate_cell_properties(_make_coo_list(masks), voxel_size)

        for col in ['label', 'area', 'x_cell', 'y_cell', 'z_cell']:
            np.testing.assert_array_equal(result[col].values, expected[col].values)

    def test_many_cells_random(self):
        """Realistic scenario: many cells spread across planes."""
        rng = np.random.default_rng(42)
        masks = np.zeros((10, 100, 100), dtype=np.uint32)

        for label in range(1, 31):
            # Each cell spans 2-4 planes
            start_plane = rng.integers(0, 7)
            n_planes = rng.integers(2, 5)
            r, c = rng.integers(5, 85, size=2)
            size = rng.integers(3, 10)
            for p in range(start_plane, min(start_plane + n_planes, 10)):
                masks[p, r:r+size, c:c+size] = label

        voxel_size = [0.28, 0.28, 0.7]
        expected = _calculate_cell_properties_regionprops(masks, voxel_size)
        result = calculate_cell_properties(_make_coo_list(masks), voxel_size)

        assert len(result) == len(expected)
        for col in ['label', 'area', 'x_cell', 'y_cell', 'z_cell']:
            np.testing.assert_array_equal(result[col].values, expected[col].values)

    def test_empty_planes_between_cells(self):
        """Empty planes between cells should not affect results."""
        masks = np.zeros((10, 20, 20), dtype=np.uint32)
        masks[0, 2:5, 2:5] = 1
        masks[1, 2:5, 2:5] = 1
        # planes 2-7 are empty
        masks[8, 10:15, 10:15] = 2
        masks[9, 10:15, 10:15] = 2

        voxel_size = [0.28, 0.28, 0.7]
        expected = _calculate_cell_properties_regionprops(masks, voxel_size)
        result = calculate_cell_properties(_make_coo_list(masks), voxel_size)

        for col in ['label', 'area', 'x_cell', 'y_cell', 'z_cell']:
            np.testing.assert_array_equal(result[col].values, expected[col].values)