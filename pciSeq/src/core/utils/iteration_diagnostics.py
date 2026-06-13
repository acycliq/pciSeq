import logging

import numpy as np

logger = logging.getLogger(__name__)


def log_iteration_diagnostics(vb, i, delta, p0, classProb_before):
    """Read-only per-iteration diagnostic logging. Changes no model state."""
    # note: this is the single largest spot-to-cell prob change, not a mean
    logger.info('Iteration %d, max spot-to-cell prob change %f' % (i, delta))
    # --- SMART LOGGING ---
    if delta > 0:
        p1 = vb.spots.parent_cell_prob
        p0_val = p0 if p0 is not None else np.zeros_like(p1)
        diffs = np.abs(p1 - p0_val)
        max_idx = np.unravel_index(np.argmax(diffs), diffs.shape)
        spot_idx = max_idx[0]
        col_idx = max_idx[1]
        gene_name = vb.spots.data.gene_name.iloc[spot_idx]

        # the cell whose assignment prob for this spot moved the most (this is what delta measures)
        moved_cell = vb.spots.parent_cell_id[spot_idx, col_idx]
        moved_prob_before = p0_val[spot_idx, col_idx]
        moved_prob_after = p1[spot_idx, col_idx]

        # the spot's single most likely parent, before and after (the winner's identity can change)
        best_col_before = np.argmax(p0_val[spot_idx])
        best_col_after = np.argmax(p1[spot_idx])
        parent_before = vb.spots.parent_cell_id[spot_idx, best_col_before]
        parent_after = vb.spots.parent_cell_id[spot_idx, best_col_after]
        prob_before = p0_val[spot_idx, best_col_before]
        prob_after = p1[spot_idx, best_col_after]

        # most likely cell type of the cell that moved the most, before and after the update
        class_prob_before = classProb_before[moved_cell]
        class_prob_after = vb.cells.classProb[moved_cell]
        class_idx_before = np.argmax(class_prob_before)
        class_idx_after = np.argmax(class_prob_after)
        class_name_before = vb.cells.class_names[class_idx_before]
        class_name_after = vb.cells.class_names[class_idx_after]
        class_p_before = class_prob_before[class_idx_before]
        class_p_after = class_prob_after[class_idx_after]

        # cell id 0 is the background "cell" (spot owned by no real cell / noise)
        def cell_label(cid):
            return "background (cell 0)" if cid == 0 else f"cell {cid}"

        logger.info(f"DIAGNOSTIC: biggest mover this iteration is spot {spot_idx} (gene {gene_name})")
        logger.info(f"DIAGNOSTIC:   most likely parent: {cell_label(parent_before)} p={prob_before:.4f} -> {cell_label(parent_after)} p={prob_after:.4f}")
        logger.info(f"DIAGNOSTIC:   largest single change: P(spot in {cell_label(moved_cell)}) went {moved_prob_before:.4f} -> {moved_prob_after:.4f} (|change|={delta:.4f})")
        logger.info(f"DIAGNOSTIC:   {cell_label(moved_cell)} most likely type: {class_name_before} (p={class_p_before:.4f}) -> {class_name_after} (p={class_p_after:.4f})")

        # ---------------------------------------------------------
        # Root-cause instrumentation for the ping-pong (read-only).
        # Reconstruct the exact softmax inputs from arrays already
        # stashed during this iteration's updates, so we can see WHICH
        # term drives the flip. No recomputation, no model change.
        if moved_cell != 0:
            # Cell-class score: wCellClass = contr (gene loglik) + log_prior + mrf
            # (matches cell_to_cellType: contr summed over genes, + prior + mrf).
            cell_contr = vb.cells.nb_contr[moved_cell].sum(axis=0)   # (nK,)
            cell_mrf = vb.cells.mrf[moved_cell]                       # (nK,)
            cell_prior = vb.cellTypes.log_prior                       # (nK,)
            cell_wtotal = cell_contr + cell_prior + cell_mrf
            a, b = np.argsort(cell_wtotal)[-2:][::-1]                    # winner, runner-up
            for k in (a, b):
                logger.info(
                    f"DIAGNOSTIC:   [class score] {cell_label(moved_cell)} {vb.cells.class_names[k]}: "
                    f"total={cell_wtotal[k]:.4f} (contr={cell_contr[k]:.4f}, prior={cell_prior[k]:.4f}, mrf={cell_mrf[k]:.4f})"
                )
            # The swing: how each component separates the top-2 classes.
            # The component whose sign flips across iterations is the driver.
            logger.info(
                f"DIAGNOSTIC:   [class swing] {vb.cells.class_names[a]} minus {vb.cells.class_names[b]}: "
                f"d_total={cell_wtotal[a] - cell_wtotal[b]:.4f} "
                f"(d_contr={cell_contr[a] - cell_contr[b]:.4f}, "
                f"d_prior={cell_prior[a] - cell_prior[b]:.4f}, "
                f"d_mrf={cell_mrf[a] - cell_mrf[b]:.4f})"
            )

            # Spot-to-cell score for the moved cell's column:
            # wSpotCell = attention(t1) + expr(t2) + cell_ineff(t3) + gene_ineff + mvn
            # (InsideCellBonus is 0 in this config). t1,t2,t3 depend on the cell's
            # class prob; mvn and gene_ineff do not. If t1/t2/t3 flip with the class
            # while mvn stays put, the spot is reacting to the class, not to geometry.
            sc = col_idx
            t1 = vb.spots.attention[spot_idx, sc]
            t2 = vb.spots.expr_fluctuations[spot_idx, sc]
            t3 = vb.spots.cell_inefficiency[spot_idx, sc]
            ge = vb.spots.gene_inefficiency[spot_idx, sc]
            mv = vb.spots.mvn_loglik_arr[spot_idx, sc]
            logger.info(
                f"DIAGNOSTIC:   [spot terms] spot {spot_idx} -> {cell_label(moved_cell)}: "
                f"sum={t1 + t2 + t3 + ge + mv:.4f} (attention={t1:.4f}, expr={t2:.4f}, "
                f"cell_ineff={t3:.4f}, gene_ineff={ge:.4f}, mvn={mv:.4f})"
            )

            # Neighbour cells that drive this cell's MRF. The MRF support is a
            # distance-weighted vote of these neighbours' classes (reconstructed
            # here exactly as in cells.calc_mrf). Watching them across
            # iterations shows which neighbours flip and produce the MRF swing.
            nbr_ids = vb.cells.nbrs['indices'][moved_cell]
            nbr_dist = vb.cells.nbrs['distances'][moved_cell]
            prox = 1.0 / nbr_dist
            prox = prox / prox.sum() * len(nbr_ids)
            nbr_cp = vb.cells.classProb[nbr_ids]
            nbr_top = np.argmax(nbr_cp, axis=1)
            nbr_desc = " | ".join(
                f"c{int(nbr_ids[j])} w={prox[j]:.2f} {vb.cells.class_names[nbr_top[j]]}({nbr_cp[j, nbr_top[j]]:.2f})"
                for j in range(len(nbr_ids))
            )
            logger.info(f"DIAGNOSTIC:   [neighbors] {cell_label(moved_cell)}: {nbr_desc}")
            # Proximity-weighted vote for the two competing classes (pre-beta, pre-A).
            # Times mrf_beta this should track the [class score] mrf values above.
            sup_a = float((prox * nbr_cp[:, a]).sum())
            sup_b = float((prox * nbr_cp[:, b]).sum())
            logger.info(
                f"DIAGNOSTIC:   [neighbor vote] {vb.cells.class_names[a]} support={sup_a:.3f} vs "
                f"{vb.cells.class_names[b]} support={sup_b:.3f} (x beta={vb.config['mrf_beta']} ~ mrf term)"
            )

            # Drill into the dominant neighbour. If a single neighbour carries most
            # of the proximity weight because it is far closer than the rest, the
            # spatial prior has collapsed to "copy that one cell". Show its raw
            # distance vs the next closest to judge whether it is a near-coincident pair.
            dom_j = int(np.argmax(prox))
            dom = int(nbr_ids[dom_j])
            next_dist = float(np.sort(nbr_dist)[1]) if len(nbr_dist) > 1 else float('nan')
            logger.info(
                f"DIAGNOSTIC:   [dominant nbr] {cell_label(moved_cell)} leans on c{dom}: "
                f"raw dist={nbr_dist[dom_j]:.3f} (next closest={next_dist:.3f}), "
                f"weight={prox[dom_j]:.2f} of {len(nbr_ids)}"
            )
            # The dominant neighbour's own neighbourhood: is moved_cell c{dom}'s
            # dominant neighbour too? (mutual-coupling test). moved_cell will appear
            # as c{moved_cell} in this list with its weight.
            dom_ids = vb.cells.nbrs['indices'][dom]
            dom_dist = vb.cells.nbrs['distances'][dom]
            dom_prox = 1.0 / dom_dist
            dom_prox = dom_prox / dom_prox.sum() * len(dom_ids)
            dom_cp = vb.cells.classProb[dom_ids]
            dom_top = np.argmax(dom_cp, axis=1)
            dom_desc = " | ".join(
                f"c{int(dom_ids[j])} d={dom_dist[j]:.2f} w={dom_prox[j]:.2f} {vb.cells.class_names[dom_top[j]]}"
                for j in range(len(dom_ids))
            )
            logger.info(f"DIAGNOSTIC:   [dominant nbr's nbrs] c{dom}: {dom_desc}")
