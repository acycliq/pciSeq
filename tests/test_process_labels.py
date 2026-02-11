"""
Tests for process_labels (sparse fastremap approach).
Verifies it produces identical results to the original dense approach.
"""

import numpy as np
import pytest
from scipy.sparse import coo_matrix
import fastremap
from pciSeq.src.preprocess.label_processing import process_labels


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