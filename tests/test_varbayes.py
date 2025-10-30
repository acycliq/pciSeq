"""
Unit tests for VarBayes algorithm core loop steps.

Tests individual methods of the VarBayes class with controlled inputs
to verify correctness of each algorithm step.
"""

import numpy as np
import pytest

from pciSeq.src.core.main import VarBayes


class TestVarBayesInitialization:
    """Test VarBayes initialization and setup."""

    def test_validate_config_with_missing_params(self):
        """Test that missing config parameters raise ValueError."""
        incomplete_config = {"max_iter": 10}  # Missing many required params

        with pytest.raises(ValueError) as excinfo:
            VarBayes._validate_config(incomplete_config)

        assert "Missing required config parameters" in str(excinfo.value)

    def test_validate_config_with_complete_params(self, base_opts):
        """Test that complete config passes validation."""
        # Should not raise any exception
        VarBayes._validate_config(base_opts)


class TestVarBayesGeneCountUpdate:
    """Test gene count update step (geneCount_upd)."""

    def test_genecount_upd_shape(self, minimal_varbayes):
        """Test that geneCount_upd produces correct output shape."""
        vb = minimal_varbayes
        vb.initialise_state()

        # Run the update
        vb.geneCount_upd()

        # Check output shape: nCells x nGenes
        assert vb.cells.geneCount.shape == (vb.nC, vb.nG)
        assert vb.cells.background_counts.shape == (vb.nG,)

    def test_genecount_upd_sum_conservation(self, minimal_varbayes):
        """Test that total gene counts are conserved."""
        vb = minimal_varbayes
        vb.initialise_state()

        # Store initial spot count
        total_spots = vb.nS

        # Run the update
        vb.geneCount_upd()

        # Sum of all gene counts should approximately equal total spots
        # (accounting for probabilistic assignments)
        total_assigned = vb.cells.geneCount.sum() + vb.cells.background_counts.sum()

        # Should be close (within 1% due to floating point)
        assert abs(total_assigned - total_spots) / total_spots < 0.01

    def test_genecount_upd_background_zero(self, minimal_varbayes):
        """Test that background (cell 0) has no gene counts."""
        vb = minimal_varbayes
        vb.initialise_state()
        vb.geneCount_upd()

        # Cell 0 (background) should have zero gene counts
        assert np.allclose(vb.cells.geneCount[0, :], 0)


class TestVarBayesCellTypeAssignment:
    """Test cell-to-cellType assignment step (cell_to_cellType)."""

    def test_cell_to_celltype_probability_sums(self, minimal_varbayes):
        """Test that cell type probabilities sum to 1 for each cell."""
        vb = minimal_varbayes

        # Initialize and run one iteration
        vb.initialise_state()
        vb.geneCount_upd()
        vb.gamma_upd()
        vb.cell_to_cellType()

        # Check that probabilities sum to 1 for each cell
        prob_sums = vb.cells.classProb.sum(axis=1)
        assert np.allclose(prob_sums, 1.0), "Cell type probabilities should sum to 1"

    def test_cell_to_celltype_output_shape(self, minimal_varbayes):
        """Test that cell_to_cellType produces correct output shape."""
        vb = minimal_varbayes
        vb.initialise_state()
        vb.geneCount_upd()
        vb.gamma_upd()
        vb.cell_to_cellType()

        # Shape should be nCells x nCellTypes
        assert vb.cells.classProb.shape == (vb.nC, vb.nK)

    def test_cell_to_celltype_probabilities_in_range(self, minimal_varbayes):
        """Test that all probabilities are in [0, 1]."""
        vb = minimal_varbayes
        vb.initialise_state()
        vb.geneCount_upd()
        vb.gamma_upd()
        vb.cell_to_cellType()

        assert np.all(vb.cells.classProb >= 0), "Probabilities should be >= 0"
        assert np.all(vb.cells.classProb <= 1), "Probabilities should be <= 1"


class TestVarBayesSpotsToCell:
    """Test spot-to-cell assignment step (spots_to_cell)."""

    def test_spots_to_cell_probability_sums(self, minimal_varbayes):
        """Test that spot-cell probabilities sum to 1 for each spot."""
        vb = minimal_varbayes
        vb.initialise_state()
        vb.geneCount_upd()
        vb.gamma_upd()
        vb.cell_to_cellType()
        vb.spots_to_cell()

        # Each spot's probabilities should sum to 1
        prob_sums = vb.spots.parent_cell_prob.sum(axis=1)
        assert np.allclose(prob_sums, 1.0), "Spot probabilities should sum to 1"

    def test_spots_to_cell_output_shape(self, minimal_varbayes):
        """Test that spots_to_cell produces correct output shape."""
        vb = minimal_varbayes
        vb.initialise_state()
        vb.geneCount_upd()
        vb.gamma_upd()
        vb.cell_to_cellType()
        vb.spots_to_cell()

        # Shape should be nSpots x nNeighbors
        assert vb.spots.parent_cell_prob.shape == (vb.nS, vb.nN)


class TestVarBayesEtaUpdate:
    """Test gene efficiency update step (eta_upd)."""

    def test_eta_upd_positive_values(self, minimal_varbayes):
        """Test that eta values remain positive."""
        vb = minimal_varbayes
        vb.initialise_state()

        # Run a few update steps
        vb.geneCount_upd()
        vb.gamma_upd()
        vb.cell_to_cellType()
        vb.spots_to_cell()
        vb.eta_upd()

        # Eta values should be positive
        assert np.all(vb.genes.eta_bar > 0), "Eta values should be positive"

    def test_eta_upd_reasonable_range(self, minimal_varbayes):
        """Test that eta values are in a reasonable range."""
        vb = minimal_varbayes
        vb.initialise_state()
        vb.geneCount_upd()
        vb.gamma_upd()
        vb.cell_to_cellType()
        vb.spots_to_cell()
        vb.eta_upd()

        # Eta values should typically be between 0.1 and 10
        # (gene efficiency shouldn't vary wildly)
        assert np.all(vb.genes.eta_bar < 100), "Eta values unusually high"
        assert np.all(vb.genes.eta_bar > 0.01), "Eta values unusually low"


class TestVarBayesConvergence:
    """Test convergence checking."""

    def test_convergence_detection(self, minimal_varbayes, base_opts):
        """Test that convergence is detected when probabilities stabilize."""
        vb = minimal_varbayes

        # Run the algorithm for a few iterations
        # Set max_iter to small value for fast testing
        vb.config["max_iter"] = 3
        vb.run()

        # Convergence tracking should populate iter_delta
        assert len(vb.iter_delta) > 0, "iter_delta should be populated"
        assert len(vb.iter_delta) <= 3, "Should not exceed max_iter"

        # Delta values should be positive (measuring change)
        assert all(delta >= 0 for delta in vb.iter_delta)


# TODO: Add tests for:
# - gamma_upd (complex delayed computations)
# - Integration tests for full loop iterations
