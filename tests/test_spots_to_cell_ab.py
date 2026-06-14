"""
Checks that spots_to_cell_numba gives the same result as spots_to_cell.

They are two implementations of the same step: spots_to_cell is the plain numpy
loop (the reference), spots_to_cell_numba is the fast numba version used by default.
This test keeps them in sync, so a change to one that diverges from the other is
caught here. Both only write parent_cell_prob and the diagnostic arrays (not their
inputs), so we can run one after the other on the same state and compare.
"""
import numpy as np


def test_spots_to_cell_numba_matches_loop(minimal_varbayes):
    vb = minimal_varbayes
    vb.initialise_state()

    # plain numpy loop, keep its result
    vb.spots_to_cell()
    p_loop = vb.spots.parent_cell_prob.copy()
    attention_loop = vb.spots.attention.copy()
    mvn_loop = vb.spots.mvn_loglik_arr.copy()

    # numba version on the same state, keep its result
    vb.spots_to_cell_numba()
    p_numba = vb.spots.parent_cell_prob.copy()
    attention_numba = vb.spots.attention.copy()
    mvn_numba = vb.spots.mvn_loglik_arr.copy()

    np.testing.assert_allclose(p_numba, p_loop, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(attention_numba, attention_loop, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(mvn_numba, mvn_loop, rtol=1e-12, atol=1e-12)