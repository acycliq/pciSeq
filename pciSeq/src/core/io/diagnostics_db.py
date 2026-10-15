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

    # --- Gather Data ---
    cells = varBayes.cells
    genes = varBayes.genes
    spots = varBayes.spots

    # format version 2 (October 2026). Version 1 kept scaled_means, a (nG, nK) float32
    # matrix per cell, 32 KB each and 800 MB on espio. It is area_factor[c] times
    # mean_expression_adj, the one matrix every cell shares, so now the matrix is
    # saved once in the metadata and the cell keeps its factor. In the same spirit
    # the spots table lost gene_inefficiency (logeta_bar per gene, ten copies of it
    # per spot; now one metadata row), the background slot of its blobs (always
    # zero) and the bonus column on a run with no InsideCellBonus, and the candidate
    # ids are an int32 blob rather than JSON text. 1.9 GB became about 0.6 GB.
    has_bonus = bool(varBayes.config.get('InsideCellBonus'))

    cursor.execute('''
        CREATE TABLE cells (
            internal_label INTEGER PRIMARY KEY,
            area_factor REAL,
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
            cov BLOB,
            scored_gene_count BLOB
        )
    ''')

    # the blobs hold nN - 1 float32, one per candidate cell; the background is the
    # last candidate and scores log_rho_bar instead, so its slot is not stored.
    # neighbor_cell_ids is the matching nN - 1 int32 internal labels.
    # z is REAL, the scaled z as the model has it. It was INTEGER until October
    # 2026 and the decimals went: z 2.875 (plane 1.15) was saved as 2, plane 0.8,
    # so the tools put one spot in ten on the plane below (see spot_z_exact).
    # x and y stay whole pixels, a fraction of a pixel changes no answer.
    cursor.execute('''
        CREATE TABLE spots (
            spot_id INTEGER PRIMARY KEY,
            gene_idx INTEGER,
            x INTEGER,
            y INTEGER,
            z REAL,
            neighbor_cell_ids BLOB,
            mvn_loglik BLOB,
            attention BLOB,
            expr_fluct BLOB,
            cell_inefficiency BLOB%s
        )
    ''' % (',\n            bonus BLOB' if has_bonus else ''))

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
    # the expected count matrix every cell shares, and each cell's factor on it:
    # scaled_exp[c] = area_factor[c] * mean_expression_adj (likelihood.scaled_exp)
    mean_expression_adj = np.asarray(varBayes.single_cell.mean_expression_adj.values, dtype=np.float32)
    area_factor = np.asarray(cells.ini_cell_props['area_factor'], dtype=np.float32)
    nC = int(cells.nC)
    nG, nK = mean_expression_adj.shape

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
        # Version 2 broke it (see the top of this function): scaled_means and
        # gene_inefficiency are gone, neighbor_cell_ids changed type. A reader of
        # version 2 reads version 1 too; a reader of version 1 refuses version 2.
        ('format_version', '2'),
        # the spots z is the full scaled z, so plane_of gives the spot's plane.
        # Files without this key have z cut to a whole number and their spot
        # planes cannot be trusted; fix_spot_z repairs them from arrow_spots
        ('spot_z_exact', '1'),
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
        # sc_mean_expression times Inefficiency, with the Zero column: what the
        # likelihood scores against. area_factor (cells table) times this row gives
        # a cell's scaled_exp back
        ('mean_expression_adj', json.dumps(mean_expression_adj.tolist())),
        # E[log eta] per gene, the gene term of every spot score. Not log(eta_bar).
        ('logeta_bar', json.dumps(np.asarray(genes.logeta_bar, dtype=np.float32).tolist())),
        # the prior at the end of the run. With cell_type_prior 'weighted', dalpha_upd
        # updates it right after the last class update, so it is one update newer than
        # the prior behind the class probabilities. On espio that is at most 0.0002 in
        # a class score and 0.00003 in a probability, too small to save the older one
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
    theta_bar_f32 = cells.theta_bar.astype(np.float32)
    gene_count_f32 = cells.geneCount.astype(np.float32)
    # the counts the last class update scored. gene_count above is one update newer:
    # spots_to_cell sets the counts again after the class update, and the Cells setter
    # keeps the ones it replaced (see Cells.geneCount). The explanation tools rebuild
    # a class score from these. None only when the counts were never updated twice
    scored = getattr(cells, '_gene_counts_prev', None)
    scored_f32 = (scored if scored is not None else cells.geneCount).astype(np.float32)
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
                float(area_factor[c]),
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
                scored_f32[c].tobytes(),
            ))
        cursor.executemany('INSERT INTO cells VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)',
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
        # the inside cell bonus the model adds before the softmax in spots_to_cell,
        # only on a run that has one
        bonus_f32 = ((spots.bonus_mask * varBayes.config['InsideCellBonus']).astype(np.float32)
                     if has_bonus else None)
        ids_i32 = np.asarray(neighbor_ids, dtype=np.int32)
        gene_idx = spots.gene_id.astype(np.int32)
        sql = ('INSERT INTO spots (spot_id, gene_idx, x, y, z, neighbor_cell_ids, mvn_loglik, '
               'attention, expr_fluct, cell_inefficiency%s) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?%s)'
               % ((', bonus', ', ?') if has_bonus else ('', '')))
        xs = spots.data['x'].astype(np.int32).to_numpy()
        ys = spots.data['y'].astype(np.int32).to_numpy()
        zs = spots.data['z'].astype(np.float32).to_numpy()

        batch_size = 10000
        for start in range(0, nS, batch_size):
            end = min(start + batch_size, nS)
            batch = []
            for i in range(start, end):
                row = [
                    int(spots.data.index[i]),
                    int(gene_idx[i]),
                    int(xs[i]), int(ys[i]), float(zs[i]),
                    ids_i32[i, :-1].tobytes(),
                    mvn_f32[i, :-1].tobytes(),
                    attn_f32[i, :-1].tobytes(),
                    expr_f32[i, :-1].tobytes(),
                    cineff_f32[i, :-1].tobytes(),
                ]
                if has_bonus:
                    row.append(bonus_f32[i, :-1].tobytes())
                batch.append(tuple(row))
            cursor.executemany(sql, batch)
            # if (end % 50000 == 0) or (end == nS):
            #     logger.info('Inserted %d/%d spots', end, nS)

    export_docs(cursor)

    conn.commit()
    conn.close()

    db_size_mb = os.path.getsize(db_path) / (1024 * 1024)
    logger.info('Saved at: %s (%.1f MB)', db_path, db_size_mb)



def shrink_diagnostics_db(src: str, dst: str) -> None:
    """Rewrite a format 1 diagnostics.db as format 2.

    TO REMOVE: this is a bridge for the format 1 files written before October
    2026. Once no run in that format is left that anyone wants to open, delete
    this function, its export in io/__init__.py, its test in
    tests/test_diagnostics_export.py (_as_format_1 and the round trip), and the
    format 1 branches of the readers: Run._cell and Run._spot_row in
    pciSeq/src/mcp/tools.py, and queryCell and querySpot in the viewer's
    electron/diagnostics.js (dbFormat === 1).

    For the runs already on disk, so they do not need a refit to get the smaller
    file: the sample run of the docs went from 1.9 GB to about 0.6 GB. Nothing is
    lost. Version 2 keeps the same numbers in another arrangement, see
    export_diagnostics: the per cell scaled_means matrix becomes area_factor times
    one shared matrix, gene_inefficiency becomes one row per gene, the background
    slot of the spot blobs goes, as does the bonus column on a run with no
    InsideCellBonus, and the candidate ids are stored as int32 rather than text.

    Parameters
    ----------
    src : str
        A format 1 diagnostics.db.
    dst : str
        Where to write the format 2 file. Must not exist.

    Raises
    ------
    ValueError
        When src is not format 1, or when its scaled_means are not area_factor
        times the shared matrix to float32 precision, which would mean the run was
        made by a pciSeq this function does not know.
    """
    if os.path.exists(dst):
        raise ValueError('%s exists already' % dst)
    old = sqlite3.connect(src)
    meta = dict(old.execute('SELECT key, value FROM metadata'))
    if meta.get('format_version') != '1':
        raise ValueError('expected a format 1 diagnostics.db, this one says %r'
                         % meta.get('format_version'))
    nG, nK, nN = int(meta['nG']), int(meta['nK']), int(meta['nN'])
    cfg = json.loads(meta['config'])
    has_bonus = bool(cfg.get('InsideCellBonus'))
    # sc_mean_expression is saved as mean_expression_adj / Inefficiency
    # (SingleCell.mean_expression), so this is the matrix the likelihood used
    m_adj = (np.asarray(json.loads(meta['sc_mean_expression']), dtype=np.float32)
             * np.float32(cfg['Inefficiency'])).astype(np.float32)
    nz = m_adj > 0

    new = sqlite3.connect(dst)
    cur = new.cursor()
    cur.execute('CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT)')

    # --- cells: the same columns with area_factor in place of scaled_means
    cols = [(r[1], r[2]) for r in old.execute('PRAGMA table_info(cells)')]
    names = [n for n, _ in cols]
    i_sm = names.index('scaled_means')
    new_cols = [('area_factor', 'REAL') if n == 'scaled_means' else (n, t) for n, t in cols]
    key = names[0]
    cur.execute('CREATE TABLE cells (%s)' % ', '.join(
        '%s %s%s' % (n, t, ' PRIMARY KEY' if n == key else '') for n, t in new_cols))
    batch = []
    for row in old.execute('SELECT * FROM cells'):
        sm = np.frombuffer(row[i_sm], dtype=np.float32).reshape(nG, nK)
        a = float(np.median(sm[nz] / m_adj[nz])) if nz.any() else 0.0
        back = (np.float32(a) * m_adj).astype(np.float32)
        if not np.allclose(back, sm, rtol=1e-5, atol=1e-6):
            raise ValueError('cell %s: scaled_means is not area_factor times the shared '
                             'matrix, this file was not made by a pciSeq this converter knows'
                             % row[0])
        batch.append(tuple(a if j == i_sm else v for j, v in enumerate(row)))
        if len(batch) == 5000:
            cur.executemany('INSERT INTO cells VALUES (%s)' % ', '.join('?' * len(cols)), batch)
            batch = []
    if batch:
        cur.executemany('INSERT INTO cells VALUES (%s)' % ', '.join('?' * len(cols)), batch)

    # --- spots: nN - 1 slots, int32 ids, no gene_inefficiency, bonus only when used
    cur.execute('''
        CREATE TABLE spots (
            spot_id INTEGER PRIMARY KEY,
            gene_idx INTEGER,
            x INTEGER,
            y INTEGER,
            z INTEGER,
            neighbor_cell_ids BLOB,
            mvn_loglik BLOB,
            attention BLOB,
            expr_fluct BLOB,
            cell_inefficiency BLOB%s
        )
    ''' % (',\n            bonus BLOB' if has_bonus else ''))
    sql = ('INSERT INTO spots VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?%s)'
           % (', ?' if has_bonus else ''))
    # gene_inefficiency is logeta_bar[gene] on every spot of the gene: one row
    logeta = np.full(nG, np.nan, dtype=np.float32)
    f32 = lambda b: np.frombuffer(b, dtype=np.float32)
    batch = []
    for (sid, g, x, y, z, ids, mvn, attn, expr, cineff, gineff, bonus) in old.execute(
            'SELECT spot_id, gene_idx, x, y, z, neighbor_cell_ids, mvn_loglik, attention, '
            'expr_fluct, cell_inefficiency, gene_inefficiency, bonus FROM spots'):
        if np.isnan(logeta[g]):
            logeta[g] = f32(gineff)[0]
        row = [sid, g, x, y, z,
               np.asarray(json.loads(ids)[:-1], dtype=np.int32).tobytes(),
               f32(mvn)[:-1].tobytes(), f32(attn)[:-1].tobytes(),
               f32(expr)[:-1].tobytes(), f32(cineff)[:-1].tobytes()]
        if has_bonus:
            row.append(f32(bonus)[:-1].tobytes())
        batch.append(tuple(row))
        if len(batch) == 20000:
            cur.executemany(sql, batch)
            batch = []
    if batch:
        cur.executemany(sql, batch)
    if np.isnan(logeta).any():
        # a gene with no spot at all never got its value; log(eta_bar) is the
        # closest thing the file holds
        eta = np.asarray(json.loads(meta['eta_bar']), dtype=np.float32)
        miss = np.isnan(logeta)
        logger.warning('%d genes have no spots, their logeta_bar is set to log(eta_bar)', miss.sum())
        logeta[miss] = np.log(eta[miss])

    # --- metadata, with the two new rows and the version
    meta['format_version'] = '2'
    meta['mean_expression_adj'] = json.dumps(m_adj.tolist())
    meta['logeta_bar'] = json.dumps(logeta.tolist())
    cur.executemany('INSERT INTO metadata VALUES (?, ?)', list(meta.items()))

    # --- docs, as they are
    has_docs = old.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='docs'").fetchone()
    cur.execute('CREATE TABLE docs (page TEXT PRIMARY KEY, text TEXT)')
    if has_docs:
        cur.executemany('INSERT INTO docs VALUES (?, ?)', old.execute('SELECT page, text FROM docs').fetchall())

    new.commit()
    new.close()
    old.close()
    logger.info('Rewrote %s (%.1f MB) as format 2 at %s (%.1f MB)', src,
                os.path.getsize(src) / 2**20, dst, os.path.getsize(dst) / 2**20)


def fix_spot_z(db_path: str, spots_dir: str = None) -> None:
    """Put the full spot z back into a diagnostics.db written before October 2026.

    TO REMOVE: a bridge like shrink_diagnostics_db. Those files saved the spot z
    as a whole number, so a spot's plane worked out from it was one plane too low
    for about one spot in ten. The viewer's arrow_spots shards of the same run
    still have the z with its decimals; this copies it over, in place, and adds
    the spot_z_exact key so the readers trust the planes again. Run it once per
    file. A format 1 file goes through shrink_diagnostics_db first.

    Parameters
    ----------
    db_path : str
        A format 2 diagnostics.db, changed in place.
    spots_dir : str, optional
        The run's arrow_spots folder. By default the one next to the diagnostics
        folder, viewer_data/arrow_spots.

    Raises
    ------
    ValueError
        When the file is not format 2, when the shards hold other spots, or when
        a shard z is not the one the file cut to a whole number, which would mean
        the shards come from another run. Nothing is written then.
    """
    import glob
    import pyarrow.feather as feather
    from ..utils.geometry import plane_of

    con = sqlite3.connect(db_path)
    try:
        meta = dict(con.execute('SELECT key, value FROM metadata'))
        if meta.get('spot_z_exact') == '1':
            logger.info('%s has the full spot z already, nothing to do', db_path)
            return
        if meta.get('format_version') != '2':
            raise ValueError('expected a format 2 diagnostics.db, this one says %r; a format '
                             '1 file goes through shrink_diagnostics_db first'
                             % meta.get('format_version'))
        if spots_dir is None:
            spots_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(db_path))),
                                     'arrow_spots')
        shards = sorted(glob.glob(os.path.join(spots_dir, 'spots_plane_*.feather')))
        if not shards:
            raise ValueError('no arrow_spots shards in %s' % spots_dir)
        tables = [feather.read_table(p, columns=['spot_id', 'z', 'plane_id']) for p in shards]
        sid = np.concatenate([t.column('spot_id').to_numpy() for t in tables]).astype(np.int64)
        z = np.concatenate([t.column('z').to_numpy() for t in tables]).astype(np.float64)
        pid = np.concatenate([t.column('plane_id').to_numpy() for t in tables]).astype(np.int64)
        order = np.argsort(sid)
        sid, z, pid = sid[order], z[order], pid[order]

        old = np.array(con.execute('SELECT spot_id, z FROM spots ORDER BY spot_id').fetchall(),
                       dtype=np.float64).reshape(-1, 2)
        if len(old) != len(sid) or not np.array_equal(old[:, 0].astype(np.int64), sid):
            raise ValueError('the shards in %s hold other spots than %s' % (spots_dir, db_path))
        # the file has int(z); the shard z is rounded to 3 decimals, so a z just
        # under a whole number can read as that number
        if not (np.abs(np.floor(z + 5e-4) - old[:, 1]) <= 1).all() or (z < old[:, 1] - 1e-3).any():
            raise ValueError('the shard z does not match the whole numbers in %s, the shards '
                             'look like another run' % db_path)
        voxel_size = json.loads(meta['config']).get('voxel_size')
        if voxel_size:
            wrong = int((plane_of(z, voxel_size) != pid).sum())
            if wrong:
                raise ValueError('%d spots would still get a plane other than their plane_id, '
                                 'nothing written' % wrong)

        with con:
            con.executemany('UPDATE spots SET z = ? WHERE spot_id = ?',
                            zip(z.tolist(), sid.tolist()))
            con.execute("INSERT INTO metadata VALUES ('spot_z_exact', '1')")
        # the update leaves the old pages behind, about 45 MB on espio; this
        # gives them back
        con.execute('VACUUM')
    finally:
        con.close()
    logger.info('Put the full spot z back into %s, %d spots', db_path, len(sid))


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
