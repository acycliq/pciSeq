"""
Checks the diagnostics database comes out with everything in it.

The dashboard and the viewer read this file, so a missing key or column is a
hard failure for them. This builds a tiny run, exports it, and looks at what
actually landed in the sqlite file rather than trusting the code reads right.
"""
import json
import os
import sqlite3

import numpy as np
import pytest

from pciSeq.src.core.io import export_diagnostics

EXPECTED_METADATA = {
    'format_version',
    'nC', 'nG', 'nK', 'rSpot', 'SpotReg', 'class_names', 'eta_bar',
    'mean_gene_reads_per_class', 'sc_mean_expression', 'mean_expression_adj',
    'logeta_bar', 'log_prior',
    'nS', 'nN', 'misread_density', 'rho_bar', 'log_rho_bar',
    'hard_misread_counts', 'hard_misread_by_plane',
    'gene_panel', 'label_map', 'gene_total_spots', 'pciSeq_provenance',
}

EXPECTED_CELL_COLUMNS = {
    'internal_label', 'area_factor', 'theta_bar', 'gene_count', 'class_prob',
    'theta', 'assigned_class_idx', 'gamma_assigned', 'mrf',
    'gene_loglik', 'runner_up_idx', 'contr_assigned', 'contr_runner_up',
    'x', 'y', 'z', 'cov',
}

# format 2: no gene_inefficiency (logeta_bar in the metadata) and bonus only on a
# run with an InsideCellBonus, which the minimal run has not
EXPECTED_SPOT_COLUMNS = {
    'spot_id', 'gene_idx', 'x', 'y', 'z', 'neighbor_cell_ids', 'mvn_loglik',
    'attention', 'expr_fluct', 'cell_inefficiency',
}


@pytest.fixture
def exported_db(minimal_varbayes, tmp_path):
    vb = minimal_varbayes
    vb.initialise_state()
    # one pass of the loop, in the order main_loop runs them, so everything the
    # export reads has actually been worked out
    vb.geneCount_upd()
    vb.rho_upd()
    vb.eta_upd()
    vb.theta_upd()
    vb.gamma_upd()
    vb.cell_to_cellType()
    vb.spots_to_cell_numba()
    export_diagnostics(vb, str(tmp_path))
    db = os.path.join(str(tmp_path), 'diagnostics', 'diagnostics.db')
    assert os.path.exists(db), 'no database was written'
    return sqlite3.connect(db)


def _columns(con, table):
    return {r[1] for r in con.execute('PRAGMA table_info(%s)' % table)}


def test_metadata_has_every_key(exported_db):
    got = {k for (k,) in exported_db.execute('select key from metadata')}
    assert EXPECTED_METADATA <= got, 'missing: %s' % (EXPECTED_METADATA - got)


def test_cells_table_columns(exported_db):
    got = _columns(exported_db, 'cells')
    assert EXPECTED_CELL_COLUMNS <= got, 'missing: %s' % (EXPECTED_CELL_COLUMNS - got)


def test_spots_table_columns(exported_db):
    got = _columns(exported_db, 'spots')
    assert EXPECTED_SPOT_COLUMNS <= got, 'missing: %s' % (EXPECTED_SPOT_COLUMNS - got)


def test_the_docs_travel_with_the_run(exported_db):
    from pciSeq.src.mcp import docs
    rows = dict(exported_db.execute('select page, text from docs'))
    assert set(rows) == set(docs.list_pages())
    assert rows['index.md'] == docs.read_page('index.md')
    # the oscillation page never leaves the repo
    assert not set(rows) & set(docs.LEFT_OUT)


def test_provenance_says_which_code_made_it(exported_db):
    (val,) = exported_db.execute(
        "select value from metadata where key='pciSeq_provenance'").fetchone()
    prov = json.loads(val)
    for k in ('version', 'branch', 'commit', 'commit_date', 'created_at'):
        assert k in prov, k


def test_the_new_blobs_are_not_empty(exported_db):
    (mrf,) = exported_db.execute('select mrf from cells limit 1').fetchone()
    assert mrf and len(mrf) > 0
    row = exported_db.execute('select neighbor_cell_ids, cell_inefficiency from spots limit 1').fetchone()
    assert row is not None and all(b is not None for b in row)


def test_format_2_holds_the_same_numbers_as_the_model(exported_db, minimal_varbayes):
    """scaled_exp is area_factor times mean_expression_adj, and the gene term of a
    spot is logeta_bar of its gene. Format 2 stores the factor, the matrix and the
    row, and this checks they give the model's own arrays back, with the spot blobs
    being the model's arrays minus the background slot."""
    vb = minimal_varbayes
    meta = dict(exported_db.execute('select key, value from metadata'))
    m_adj = np.asarray(json.loads(meta['mean_expression_adj']), dtype=np.float32)
    logeta = np.asarray(json.loads(meta['logeta_bar']), dtype=np.float32)
    np.testing.assert_allclose(m_adj, vb.single_cell.mean_expression_adj.values, rtol=1e-6)
    np.testing.assert_allclose(logeta, vb.genes.logeta_bar, rtol=1e-6)
    for c, a in exported_db.execute('select internal_label, area_factor from cells'):
        np.testing.assert_allclose(np.float32(a) * m_adj, vb.scaled_exp[c], rtol=1e-5)
    f32 = lambda b: np.frombuffer(b, dtype=np.float32)
    nN = int(meta['nN'])
    for i, (sid, g, ids, mvn, cineff) in enumerate(exported_db.execute(
            'select spot_id, gene_idx, neighbor_cell_ids, mvn_loglik, cell_inefficiency from spots')):
        ids = np.frombuffer(ids, dtype=np.int32)
        assert len(ids) == nN - 1 and len(f32(mvn)) == nN - 1
        np.testing.assert_array_equal(ids, vb.spots.parent_cell_id[i, :-1])
        np.testing.assert_allclose(f32(mvn), vb.spots.mvn_loglik_arr[i, :-1], rtol=1e-6)
        np.testing.assert_allclose(f32(cineff), vb.spots.cell_inefficiency[i, :-1], rtol=1e-6)
        assert logeta[g] == pytest.approx(vb.spots.gene_inefficiency[i, 0], rel=1e-6)


def _as_format_1(con2, path):
    """A format 1 copy of a format 2 database, laid out as export_diagnostics wrote it
    before October 2026: scaled_means per cell, gene_inefficiency and bonus per spot,
    the background slot in every blob, the ids as JSON text."""
    meta = dict(con2.execute('select key, value from metadata'))
    nG, nK, nN = int(meta['nG']), int(meta['nK']), int(meta['nN'])
    m_adj = np.asarray(json.loads(meta['mean_expression_adj']), dtype=np.float32)
    logeta = np.asarray(json.loads(meta['logeta_bar']), dtype=np.float32)
    con1 = sqlite3.connect(path)
    con1.execute('create table metadata (key text primary key, value text)')
    meta['format_version'] = '1'
    del meta['mean_expression_adj'], meta['logeta_bar']
    con1.executemany('insert into metadata values (?, ?)', list(meta.items()))
    cols = [r[1] for r in con2.execute('pragma table_info(cells)')]
    con1.execute('create table cells (%s)' % ', '.join(
        'internal_label integer primary key' if c == 'internal_label'
        else 'scaled_means blob' if c == 'area_factor' else c + ' blob' for c in cols))
    rows = []
    for row in con2.execute('select * from cells'):
        row = list(row)
        row[cols.index('area_factor')] = (np.float32(row[cols.index('area_factor')]) * m_adj).tobytes()
        rows.append(tuple(row))
    con1.executemany('insert into cells values (%s)' % ', '.join('?' * len(cols)), rows)
    con1.execute('create table spots (spot_id integer primary key, gene_idx integer, x integer, '
                 'y integer, z integer, neighbor_cell_ids text, mvn_loglik blob, attention blob, '
                 'expr_fluct blob, cell_inefficiency blob, gene_inefficiency blob, bonus blob)')
    z = np.zeros(1, dtype=np.float32)
    grow = lambda b: np.append(np.frombuffer(b, dtype=np.float32), z).tobytes()
    rows = []
    for sid, g, x, y, zz, ids, mvn, attn, expr, cineff in con2.execute(
            'select spot_id, gene_idx, x, y, z, neighbor_cell_ids, mvn_loglik, attention, '
            'expr_fluct, cell_inefficiency from spots'):
        ids = np.frombuffer(ids, dtype=np.int32).tolist() + [0]
        gi = np.append(np.full(nN - 1, logeta[g], dtype=np.float32), z).tobytes()
        rows.append((sid, g, x, y, zz, json.dumps(ids), grow(mvn), grow(attn), grow(expr),
                     grow(cineff), gi, np.zeros(nN, dtype=np.float32).tobytes()))
    con1.executemany('insert into spots values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', rows)
    con1.execute('create table docs (page text primary key, text text)')
    con1.executemany('insert into docs values (?, ?)', con2.execute('select page, text from docs').fetchall())
    con1.commit()
    con1.close()


def test_the_converter_and_both_formats_read_the_same(exported_db, tmp_path):
    """A format 1 file, converted, is byte for byte the format 2 file in every
    number, and the reader gives the same answers from either."""
    from pciSeq.src.core.io import shrink_diagnostics_db
    from pciSeq.src.mcp.tools import open_run
    v1_dir = tmp_path / 'v1' / 'diagnostics'; v1_dir.mkdir(parents=True)
    _as_format_1(exported_db, str(v1_dir / 'diagnostics.db'))
    v2_dir = tmp_path / 'v2' / 'diagnostics'; v2_dir.mkdir(parents=True)
    shrink_diagnostics_db(str(v1_dir / 'diagnostics.db'), str(v2_dir / 'diagnostics.db'))
    con2 = sqlite3.connect(v2_dir / 'diagnostics.db')
    assert dict(con2.execute("select key, value from metadata"))['format_version'] == '2'
    assert {r[1] for r in con2.execute('pragma table_info(spots)')} == EXPECTED_SPOT_COLUMNS
    for table, key in (('cells', 'internal_label'), ('spots', 'spot_id')):
        a = exported_db.execute('select * from %s order by %s' % (table, key)).fetchall()
        b = con2.execute('select * from %s order by %s' % (table, key)).fetchall()
        assert len(a) == len(b) and len(a) > 0
        for ra, rb in zip(a, b):
            for va, vb in zip(ra, rb):
                if isinstance(va, float):
                    assert va == pytest.approx(vb, rel=1e-6)
                else:
                    assert va == vb
    # the same questions, the same answers, whichever file
    r1, r2 = open_run(tmp_path / 'v1'), open_run(tmp_path / 'v2')
    assert r1.format_version == 1 and r2.format_version == 2
    cells = [c for (c,) in con2.execute('select internal_label from cells where internal_label > 0 limit 3')]
    for c in cells:
        lab = r2.to_external(c)
        e1, e2 = r1.explain_cell(lab), r2.explain_cell(lab)
        assert e1['genes_favouring_assigned'] == e2['genes_favouring_assigned']
        assert e1['genes_favouring_compared'] == e2['genes_favouring_compared']
    for (sid,) in con2.execute('select spot_id from spots limit 5'):
        assert r1.explain_spot(sid)['candidates'] == r2.explain_spot(sid)['candidates']


def test_the_saved_class_score_gives_back_the_stored_probabilities(exported_db):
    """gene_loglik is the run's own gene term, from its last class update, so with
    the prior and the spatial term it has to give the stored class_prob back through
    a softmax, for every cell. That is the point of saving it: an explanation read
    from it is the run's, whatever pciSeq reads the file years later."""
    from scipy.special import softmax
    meta = dict(exported_db.execute('select key, value from metadata'))
    nG, nK = int(meta['nG']), int(meta['nK'])
    log_prior = np.asarray(json.loads(meta['log_prior']), dtype=np.float32)
    rows = exported_db.execute(
        'select class_prob, mrf, gene_loglik, assigned_class_idx, runner_up_idx, '
        'contr_assigned, contr_runner_up from cells').fetchall()
    assert rows
    f32 = lambda b: np.frombuffer(b, dtype=np.float32)
    for cp, mrf, gl, a, r, ca, cr in rows:
        cp, mrf, gl, ca, cr = f32(cp), f32(mrf), f32(gl), f32(ca), f32(cr)
        assert gl.shape == (nK,) and ca.shape == (nG,) and cr.shape == (nG,)
        np.testing.assert_allclose(softmax(gl + log_prior + mrf), cp, atol=1e-5)
        # the per gene terms add up to the saved totals
        assert ca.sum() == pytest.approx(gl[a], rel=1e-4, abs=1e-3)
        if r >= 0:
            assert r != a and cp[r] > 0 and cp[r] == np.sort(cp)[-2]
            assert cr.sum() == pytest.approx(gl[r], rel=1e-4, abs=1e-3)
        else:
            assert (np.delete(cp, a) == 0).all()


def test_the_format_version_is_2_and_a_future_file_is_refused(exported_db, tmp_path):
    """A reader can trust any file at or below its own version and must refuse one
    from above, rather than guess."""
    (v,) = exported_db.execute("select value from metadata where key='format_version'").fetchone()
    assert v == '2'
    from pciSeq.src.mcp.tools import open_run
    run = open_run(tmp_path)
    assert run.summary()['format_version'] == 2
    assert run.run_info()['format_version'] == 2
    # a file written by a future pciSeq that had to break the promise again
    import sqlite3 as sq
    db = tmp_path / 'diagnostics' / 'diagnostics.db'
    con = sq.connect(db)
    con.execute("update metadata set value='3' where key='format_version'")
    con.commit(); con.close()
    with pytest.raises(NotImplementedError, match='format version 3'):
        open_run(tmp_path)


def test_the_centroids_land_in_the_cells_table(exported_db, minimal_varbayes):
    """x, y and the scaled z per cell, so a reader places cells from this file
    alone, without the viewer's arrow files."""
    want = minimal_varbayes.cells.centroid
    got = exported_db.execute(
        'select internal_label, x, y, z from cells order by internal_label').fetchall()
    assert len(got) == len(want)
    for i, x, y, z in got[1:]:          # row 0 is the background dummy
        assert x == pytest.approx(want['x'].iloc[i]) and y == pytest.approx(want['y'].iloc[i])
        assert z == pytest.approx(want['z'].iloc[i])
