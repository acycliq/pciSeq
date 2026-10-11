"""The sqlite database the viewer reads for its per cell and per spot panels."""

import os
import json
import sqlite3
import logging
import numpy as np
from typing import Any

logger = logging.getLogger(__name__)



def export_diagnostics(varBayes: Any, output_dir: str) -> None:
    """Export diagnostics data (check_cell and check_spot) to a single SQLite database.

    Writes to {output_dir}/diagnostics/diagnostics.db

    Tables:
      - metadata: key-value pairs (including JSON arrays)
      - cells: per-cell diagnostic data
      - spots: per-spot diagnostic data
      - docs: the documentation pages, one row per page (see export_docs)

    Args:
        varBayes: Fitted VarBayes object
        output_dir: Base data directory
    """
    diagnostics_dir = os.path.join(output_dir, 'diagnostics')
    os.makedirs(diagnostics_dir, exist_ok=True)

    db_path = os.path.join(diagnostics_dir, 'diagnostics.db')
    if os.path.exists(db_path):
        os.remove(db_path)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # --- Create Tables ---
    cursor.execute('CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT)')

    cursor.execute('''
        CREATE TABLE cells (
            internal_label INTEGER PRIMARY KEY,
            scaled_means BLOB,
            theta_bar BLOB,
            gene_count BLOB,
            class_prob BLOB,
            theta REAL,
            assigned_class_idx INTEGER,
            gamma_assigned BLOB,
            mrf BLOB,
            neighbours BLOB,
            x REAL,
            y REAL,
            z REAL,
            gene_loglik BLOB,
            runner_up_idx INTEGER,
            contr_assigned BLOB,
            contr_runner_up BLOB,
            cov BLOB
        )
    ''')

    cursor.execute('''
        CREATE TABLE spots (
            spot_id INTEGER PRIMARY KEY,
            gene_idx INTEGER,
            x INTEGER,
            y INTEGER,
            z INTEGER,
            neighbor_cell_ids TEXT,
            mvn_loglik BLOB,
            attention BLOB,
            expr_fluct BLOB,
            cell_inefficiency BLOB,
            gene_inefficiency BLOB,
            bonus BLOB
        )
    ''')

    # --- Gather Data ---
    cells = varBayes.cells
    genes = varBayes.genes
    spots = varBayes.spots

    # Label Map
    label_map = {}
    if varBayes.config.get('label_map'):
        label_map = {str(k): int(v) for k, v in varBayes.config['label_map'].items()}

    # Gene Panel
    gene_panel = genes.gene_panel.tolist()

    # Misread Density
    misread_series = genes.misread_density
    if hasattr(misread_series, 'to_dict'):
        misread_dict = {str(k): float(v) for k, v in misread_series.to_dict().items()}
    else:
        logger.error("Diagnostics export skipped: 'misread_density' is missing or invalid.")
        return

    # The misread density the model actually learned, as opposed to the prior above.
    rho_bar = genes.rho_bar
    rho_bar_dict = ({str(g): float(v) for g, v in zip(genes.gene_panel, rho_bar)}
                    if rho_bar is not None else {})

    # E[log rho] per gene. This is the term the background column in spots_to_cell
    # uses, not log(rho_bar), because E[log rho] is not log(E[rho]).
    log_rho_bar = genes.log_rho_bar
    log_rho_bar_dict = ({str(g): float(v) for g, v in zip(genes.gene_panel, log_rho_bar)}
                        if log_rho_bar is not None else {})

    # Spots whose best guess is the background column, counted per gene and per
    # gene per plane.
    hard_misread_counts = []
    hard_misread_by_plane = {}
    try:
        prob = spots.parent_cell_prob
        if prob is not None and len(prob):
            is_misread = np.argmax(prob, axis=1) == (prob.shape[1] - 1)
            hard_misread_counts = np.bincount(
                spots.gene_id, weights=is_misread.astype(np.float32),
                minlength=len(genes.gene_panel)).astype(int).tolist()
            if 'plane_id' in spots.data.columns:
                plane_ids = spots.data['plane_id'].to_numpy()
                for idx in np.where(is_misread)[0]:
                    gene = genes.gene_panel[spots.gene_id[idx]]
                    plane = int(plane_ids[idx])
                    hard_misread_by_plane.setdefault(gene, {})
                    hard_misread_by_plane[gene][plane] = hard_misread_by_plane[gene].get(plane, 0) + 1
    except Exception as e:
        logger.warning('Could not compute hard misread counts: %s', e)

    # --- Populate Metadata ---
    # Compute scaled_means for metadata nC (and for cells table)
    # logger.info('Computing scaled_exp for diagnostics export...')
    scaled_means = varBayes.scaled_exp
    nC, nG, nK = scaled_means.shape

    nS = spots.nS
    # Check if we can get nN (needs neighbor_ids)
    neighbor_ids = spots.parent_cell_id
    nN = 0
    if neighbor_ids is not None:
        nN = neighbor_ids.shape[1]

    meta_items = [
        # The format of this file. A run is often opened months or years after the
        # fit, so from version 1 on the promise is: columns and keys are only ever
        # ADDED, never renamed, removed or given a new meaning. A reader written for
        # version 1 then reads every later file. Bump the number only when that
        # promise has to be broken, and expect readers to refuse what they do not
        # know. Runs from before this key are pre-1 and unsupported: regenerate.
        ('format_version', '1'),
        # Cell-related
        ('nC', str(nC)),
        ('nG', str(nG)),
        ('nK', str(nK)),
        ('rSpot', str(float(varBayes.config['rSpot']))),
        ('SpotReg', str(float(varBayes.config['SpotReg']))),
        # the mean cell radius, in pixels of the xy plane. It sets the Gaussian of
        # every cell (cov = mcr^2 I, the Gaussian update is not run) and it was
        # only recoverable from cellData.sphere_scale / 3, which nobody would guess
        ('mcr', str(float(cells.mcr))),
        ('class_names', json.dumps(cells.class_names.tolist())),
        ('eta_bar', json.dumps(genes.eta_bar.astype(np.float32).tolist())),
        ('mean_gene_reads_per_class', json.dumps(cells.mean_gene_reads_per_class().astype(np.float32).tolist())),
        ('sc_mean_expression', json.dumps(varBayes.single_cell.mean_expression.values.astype(np.float32).tolist())),
        ('log_prior', json.dumps(varBayes.cellTypes.log_prior.astype(np.float32).tolist())),

        # Spot-related
        ('nS', str(int(nS))),
        ('nN', str(int(nN))),
        ('misread_density', json.dumps(misread_dict)),
        ('rho_bar', json.dumps(rho_bar_dict)),
        ('log_rho_bar', json.dumps(log_rho_bar_dict)),
        ('hard_misread_counts', json.dumps(hard_misread_counts)),
        ('hard_misread_by_plane', json.dumps(hard_misread_by_plane)),

        # Shared
        ('gene_panel', json.dumps(gene_panel)),
        ('label_map', json.dumps(label_map)),

        # The resolved settings and how the run ended, so a reader of this file can
        # say what produced it without the pickle. label_map has its own row above.
        # voxel_size rides inside config, which is what turns the scaled z of the
        # spots table back into a plane number.
        ('config', json.dumps({k: v for k, v in varBayes.config.items() if k != 'label_map'},
                              default=str)),
        ('run', json.dumps({
            'iterations': int(varBayes.iter_num) + 1 if varBayes.iter_num is not None else None,
            'converged': bool(varBayes.has_converged),
            'delta': [float(d) for d in varBayes.iter_delta],
        })),
    ]

    # Per-gene observed spot counts (for η scatter in dashboard)
    try:
        if hasattr(spots, 'counts_per_gene') and spots.counts_per_gene is not None:
            gene_total_spots = spots.counts_per_gene.astype(int).tolist()
        else:
            vc = spots.data['gene_name'].value_counts()
            gene_total_spots = [int(vc.get(g, 0)) for g in gene_panel]
        meta_items.append(('gene_total_spots', json.dumps(gene_total_spots)))
    except Exception:
        pass

    # which pciSeq made this. version, branch, commit, commit_date, created_at,
    # serialised_at, hostname, os, python version, package versions.
    meta_items.append(('pciSeq_provenance', json.dumps(getattr(varBayes, 'metadata', {}))))

    cursor.executemany('INSERT INTO metadata VALUES (?, ?)', meta_items)
    # logger.info('Inserted %d metadata entries', len(meta_items))

    # --- Populate Cells Table ---
    scaled_means_f32 = scaled_means.astype(np.float32)
    theta_bar_f32 = cells.theta_bar.astype(np.float32)
    gene_count_f32 = cells.geneCount.astype(np.float32)
    class_prob_f32 = cells.classProb.astype(np.float32)

    # Dashboard columns: theta, assigned_class_idx, gamma_assigned
    # theta = sum_k zeta[c,k] * theta_bar[c,k]
    theta_scalar = np.einsum('ck,ck->c', class_prob_f32, theta_bar_f32).astype(np.float32)
    assigned_class_idx = np.argmax(class_prob_f32, axis=1).astype(np.int32)

    # gamma_assigned[c, :] = gamma_bar[c, :, assigned_class[c]]
    gamma_bar = varBayes.spots.gamma_bar.astype(np.float32)  # (nC, nG, nK)
    gamma_assigned = gamma_bar[np.arange(nC), :, assigned_class_idx]   # (nC, nG)

    # the mrf term the last class update used. None if the run never got that
    # far, write zeros then so the column is still the right shape.
    mrf_f32 = (cells.mrf.astype(np.float32) if cells.mrf is not None
               else np.zeros((nC, nK), dtype=np.float32))

    # the covariance of each cell's Gaussian, d x d row by row, d = 2 or 3. With the
    # cell geometry update it is fitted per cell, and before this column the only
    # trace of it was cellData's sphere_scale and sphere_rotation, 3D only and in
    # Euler angles. mcr^2 I for every cell when the update is not run.
    cov_f32 = np.ascontiguousarray(cells.cov, dtype=np.float32)

    # the cells the mrf term listens to, nNeighbors per cell, nearest first, the cell
    # itself left out. Internal labels like the key column, so a reader has to go through
    # label_map to get segmentation labels. Only in the pickle before, which is far
    # too big to open just for these.
    nbrs = cells.nbrs['indices'] if cells.nbrs is not None else np.zeros((nC, 0))
    nbrs_i32 = np.ascontiguousarray(nbrs, dtype=np.int32)

    # the centroid, in image pixels for x and y and the anisotropy scaled z the
    # model works in, so a tool can place a cell (and its plane) from this file
    # alone. Row 0 is the background dummy at -sys.maxsize, stored as it is.
    cx = cells.centroid['x'].to_numpy(dtype=np.float64)
    cy = cells.centroid['y'].to_numpy(dtype=np.float64)
    cz = cells.centroid['z'].to_numpy(dtype=np.float64)

    # The class score as the run itself computed it, from the last class update
    # (cells.nb_contr, the per cell, gene and class negative binomial terms). A run
    # is often opened months or years later by a newer pciSeq; saved here, the
    # explanation reads the run's own numbers instead of recomputing them with
    # formulas that may have changed since. gene_loglik + log_prior + mrf gives the
    # stored class_prob back through a softmax.
    #   gene_loglik      per cell and class, the sum over genes
    #   runner_up_idx    the second most probable class, -1 when no other class
    #                    has any probability at all
    #   contr_assigned   per gene, for the assigned class
    #   contr_runner_up  per gene, for the runner up (zeros when there is none)
    # The per gene terms of the other classes are left out, all of them would be
    # nC x nG x nK, about 800 MB on espio.
    nb_contr = getattr(cells, 'nb_contr', None)
    rows = np.arange(nC)
    order = np.argsort(-class_prob_f32, axis=1)
    second = np.where(order[:, 0] == assigned_class_idx, order[:, 1], order[:, 0]) if nK > 1 \
        else np.full(nC, -1)
    has_second = (second >= 0) & (class_prob_f32[rows, np.maximum(second, 0)] > 0)
    runner_up_idx = np.where(has_second, second, -1).astype(np.int32)
    if nb_contr is not None:
        gene_loglik_f32 = nb_contr.sum(axis=1).astype(np.float32)                  # (nC, nK)
        contr_assigned = nb_contr[rows, :, assigned_class_idx].astype(np.float32)   # (nC, nG)
        contr_runner_up = np.where(has_second[:, None],
                                   nb_contr[rows, :, np.maximum(second, 0)], 0).astype(np.float32)
    else:
        # the run never got to a class update, nothing to save
        gene_loglik_f32 = np.zeros((nC, nK), dtype=np.float32)
        contr_assigned = contr_runner_up = np.zeros((nC, nG), dtype=np.float32)

    batch_size = 10000
    for batch_start in range(0, nC, batch_size):
        batch_end = min(batch_start + batch_size, nC)
        batch_data = []
        for c in range(batch_start, batch_end):
            batch_data.append((
                c,
                scaled_means_f32[c].tobytes(),
                theta_bar_f32[c].tobytes(),
                gene_count_f32[c].tobytes(),
                class_prob_f32[c].tobytes(),
                float(theta_scalar[c]),
                int(assigned_class_idx[c]),
                gamma_assigned[c].tobytes(),
                mrf_f32[c].tobytes(),
                nbrs_i32[c].tobytes(),
                float(cx[c]), float(cy[c]), float(cz[c]),
                gene_loglik_f32[c].tobytes(),
                int(runner_up_idx[c]),
                contr_assigned[c].tobytes(),
                contr_runner_up[c].tobytes(),
                cov_f32[c].tobytes(),
            ))
        cursor.executemany('INSERT INTO cells VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
                           batch_data)
        # if (batch_end % 10000 == 0) or (batch_end == nC):
        #     logger.info('Inserted %d/%d cells', batch_end, nC)

    # --- Populate Spots Table ---
    if (spots.mvn_loglik_arr is None or spots.attention is None
            or spots.expr_fluctuations is None or spots.cell_inefficiency is None
            or spots.gene_inefficiency is None or neighbor_ids is None):
        logger.warning('check_spot data missing; spots table will be empty.')
    else:
        mvn_f32 = spots.mvn_loglik_arr.astype(np.float32)
        attn_f32 = spots.attention.astype(np.float32)
        expr_f32 = spots.expr_fluctuations.astype(np.float32)
        cineff_f32 = spots.cell_inefficiency.astype(np.float32)
        gineff_f32 = spots.gene_inefficiency.astype(np.float32)
        # the inside cell bonus the model adds before the softmax in spots_to_cell
        bonus_f32 = (spots.bonus_mask * varBayes.config['InsideCellBonus']).astype(np.float32)
        gene_idx = spots.gene_id.astype(np.int32)
        xs = spots.data['x'].astype(np.int32).to_numpy()
        ys = spots.data['y'].astype(np.int32).to_numpy()
        zs = spots.data['z'].astype(np.int32).to_numpy()

        batch_size = 10000
        for start in range(0, nS, batch_size):
            end = min(start + batch_size, nS)
            batch = []
            for i in range(start, end):
                neigh_json = json.dumps(list(map(int, neighbor_ids[i].tolist())))
                batch.append((
                    int(spots.data.index[i]),
                    int(gene_idx[i]),
                    int(xs[i]), int(ys[i]), int(zs[i]),
                    neigh_json,
                    mvn_f32[i].tobytes(),
                    attn_f32[i].tobytes(),
                    expr_f32[i].tobytes(),
                    cineff_f32[i].tobytes(),
                    gineff_f32[i].tobytes(),
                    bonus_f32[i].tobytes(),
                ))
            cursor.executemany('''
                INSERT INTO spots (spot_id, gene_idx, x, y, z, neighbor_cell_ids, mvn_loglik, attention, expr_fluct, cell_inefficiency, gene_inefficiency, bonus)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''', batch)
            # if (end % 50000 == 0) or (end == nS):
            #     logger.info('Inserted %d/%d spots', end, nS)

    export_docs(cursor)

    conn.commit()
    conn.close()

    db_size_mb = os.path.getsize(db_path) / (1024 * 1024)
    logger.info('Saved at: %s (%.1f MB)', db_path, db_size_mb)



def export_docs(cursor: Any) -> None:
    """Save the documentation pages (markdown only) into the docs table.

    So the docs travel with the run: the viewer's chat reads the pages of the pciSeq
    that made the numbers, with no internet and no drift. About 280 kB. The pages
    come from docs.list_pages(), the repo checkout or the copy packed into a pip
    install, and pages kept out on purpose (docs.LEFT_OUT) stay out. When no docs
    are found the table is still made, empty, and the agent says the documentation
    is not available.
    """
    from pciSeq.src.mcp import docs

    cursor.execute('CREATE TABLE docs (page TEXT PRIMARY KEY, text TEXT)')
    pages = docs.list_pages()
    if not pages:
        logger.warning('No documentation pages found, the docs table is left empty')
        return
    cursor.executemany('INSERT INTO docs VALUES (?, ?)',
                       [(p, docs.read_page(p)) for p in pages])


def export_db_tables(out_dir: str, con: Any) -> None:
    """Export all database tables to CSV files.

    Args:
        out_dir: Output directory for CSV files
        con: Database connection object
    """
    tables = con.get_db_tables()
    for table in tables:
        export_db_table(table, out_dir, con)



def export_db_table(table_name: str, out_dir: str, con: Any) -> None:
    """Export single database table to CSV.

    Args:
        table_name: Name of table to export
        out_dir: Output directory
        con: Database connection object
    """
    df = con.from_redis(table_name)
    fname = os.path.join(out_dir, table_name + '.csv')
    df.to_csv(fname, index=False)
    logger.info('Saved at %s', fname)
