"""The agent tools, checked against the run they read.

Every tool here is cross-checked against something the run wrote independently:
explain_cell against the stored class probabilities, explain_spot against
geneData.neighbour, cell_counts against cellData.CellGeneCount, spots_in_cell
against geneData.inside_cell. If a tool drifts from the model it fails here.

Uses the 16 cell fixture from test_label_identifiers, labels 101..116, so the
label translation is exercised too.
"""
import numpy as np
import pandas as pd
import pytest
from scipy.special import softmax

from pciSeq.src.mcp.tools import open_run, Run
from tests.test_label_identifiers import _run

LABELS = list(range(101, 117))


@pytest.fixture(scope='module')
def fitted(tmp_path_factory):
    """One run, saved, plus the frames it returned.

    Run long enough to converge: explain_cell recomputes the class score with the
    final eta, and eta is updated after the last class update, so the recomputation
    only matches the stored probabilities exactly once the run has settled. At two
    iterations the error is 1e-2, at ten it is 1e-7."""
    tmp = tmp_path_factory.mktemp('mcp')
    rng = np.random.default_rng(7)
    # tol 0.5 would call it converged after two iterations whatever max_iter says
    cellData, geneData = _run(rng, tmp, save=True, max_iter=40, tol=0.001)
    return open_run(tmp), cellData.set_index('Cell_Num'), geneData.set_index('spot_id')


# ------------------------------------------------------------ opening it

def test_open_run_finds_the_db_from_the_run_folder(fitted):
    run, _, _ = fitted
    assert run.path.name == 'diagnostics.db'


def test_summary_does_not_count_the_background(fitted):
    run, cellData, _ = fitted
    s = run.summary()
    assert s['cells'] == 16 == len(cellData)
    assert s['labels_were_renumbered'] is True
    assert s['has_containment'] is True


def test_labels_round_trip(fitted):
    run, _, _ = fitted
    for lab in LABELS:
        assert run.to_external(run.to_internal(lab)) == lab


# ---------------------------------------------------------- explain_cell

def test_explain_cell_reproduces_the_stored_class_probabilities(fitted):
    """The whole point. The recomputed score must give back what the model stored,
    otherwise the gene by gene story is not the model's story."""
    run, _, _ = fitted
    log_prior = np.asarray(run._meta('log_prior', parse=True), dtype=np.float32)
    for lab in LABELS:
        c = run._cell(lab)
        mine = softmax(run._nb_contributions(c).sum(axis=0) + log_prior + c['mrf'])
        assert np.allclose(mine, c['class_prob'], atol=1e-4), lab


def test_the_rebuild_is_the_models_own_score_gene_by_gene(fitted):
    """The explanation tools rebuild the class score of any class from the saved
    inputs (_nb_contributions here, nbContributions in the viewer). It must be the
    model's own: gene by gene for the classes whose per gene numbers the run saved,
    and in total for every class. If cell_to_cellType or the likelihood it calls
    changes and the rebuild is not changed with it, this fails."""
    run, _, _ = fitted
    assert run._has_scored_counts, 'the run did not save the counts its class update scored'
    for lab in LABELS:
        c = run._cell(lab)
        contr = run._nb_contributions(c).astype(np.float64)
        saved = run._saved_score(c['row'])
        a = c['assigned_class_idx']
        np.testing.assert_allclose(contr[:, a], saved['contr_assigned'], rtol=1e-4, atol=1e-3)
        if saved['runner_up_idx'] >= 0:
            np.testing.assert_allclose(contr[:, saved['runner_up_idx']], saved['contr_runner_up'],
                                       rtol=1e-4, atol=1e-3)
        np.testing.assert_allclose(contr.sum(axis=0), saved['gene_loglik'], rtol=1e-4, atol=0.05)
        # and so explain_cell gives the genes against every class, none refused
        for k in run.class_names:
            if k != run.class_names[a]:
                assert 'detail_note' not in run.explain_cell(lab, k), (lab, k)


def test_explain_cell_agrees_with_cellData(fitted):
    run, cellData, _ = fitted
    for lab in LABELS:
        e = run.explain_cell(lab)
        assert e['assigned'] == cellData.loc[lab, 'ClassName'][0]
        # cellData cuts Prob to 3 decimals (floor, not round), so the two can be a
        # full 0.001 apart, plus float noise on top. 0.9999 vs 0.999 is fine
        assert e['prob_assigned'] == pytest.approx(cellData.loc[lab, 'Prob'][0], abs=1.01e-3)
        assert e['assigned'] != e['compared_with']


def test_shared_genes_are_big_counts_that_do_not_separate_the_classes(fitted):
    """The finalists share them: among the cell's ten biggest counts, each fitting
    both classes within a factor of e, biggest first, at most five."""
    run, _, _ = fitted
    for lab in LABELS:
        e = run.explain_cell(lab)
        shared = e['shared_genes']
        assert len(shared) <= 5
        assert all(abs(g['diff']) < 1.0 for g in shared)
        counts = [g['counts'] for g in shared]
        assert counts == sorted(counts, reverse=True)
        c = run._cell(lab)
        ten_biggest = set(run.gene_panel[np.argsort(-c['gene_count'])[:10]])
        assert {g['gene'] for g in shared} <= ten_biggest
        assert 'finalists' in e['shared_genes_are']


def test_no_runner_up_asks_for_a_class_rather_than_picking_one(fitted, monkeypatch):
    """A cell whose class holds all the probability has no runner up. argsort would
    return the first class of the list and the story would call it the closest
    alternative; the tool asks instead. Seen on espio cells 2632 and 24454."""
    run, _, _ = fitted
    real = run._cell

    def one_hot(label):
        c = dict(real(label))
        p = np.zeros_like(c['class_prob'])
        p[c['assigned_class_idx']] = 1.0
        c['class_prob'] = p
        return c

    monkeypatch.setattr(run, '_cell', one_hot)
    with pytest.raises(ValueError, match='no runner up'):
        run.explain_cell(101)
    # naming a class still works
    other = [k for k in run.class_names if k != run.cell(101)['classes'][0]['class']][0]
    assert run.explain_cell(101, vs_class=other)['compared_with'] == other


def test_explain_cell_against_a_named_class(fitted):
    run, _, _ = fitted
    e = run.explain_cell(101)
    other = 'TypeA' if e['assigned'] != 'TypeA' else 'TypeB'
    assert run.explain_cell(101, vs_class=other)['compared_with'] == other
    with pytest.raises(ValueError):
        run.explain_cell(101, vs_class=e['assigned'])
    with pytest.raises(ValueError):
        run.explain_cell(101, vs_class='NoSuchClass')


# ---------------------------------------------------------- explain_spot

def test_explain_spot_agrees_with_geneData(fitted):
    """The softmax of the stored terms has to land on the same cell geneData did,
    and the candidates have to come out as segmentation labels."""
    run, _, geneData = fitted
    for sid in geneData.index[:200]:
        e = run.explain_spot(sid)
        probs = [c['prob'] for c in e['candidates']]
        assert sum(probs) == pytest.approx(1.0, abs=1e-5)

        expected = geneData.loc[sid, 'neighbour']
        assert e['assigned_to'] == ('background' if expected == 0 else expected)

        cells = [c['cell'] for c in e['candidates'] if c['cell'] != 'background']
        assert set(cells) <= set(LABELS)


def test_explain_spot_gene_matches(fitted):
    run, _, geneData = fitted
    sid = geneData.index[0]
    assert run.explain_spot(sid)['gene'] == geneData.loc[sid, 'gene_name']


# ----------------------------------------------------------- cell_counts

def test_cell_counts_matches_cellData_and_says_it_is_soft(fitted):
    run, cellData, _ = fitted
    for lab in LABELS:
        cc = run.cell_counts(lab)
        assert cc['total_counts'] == pytest.approx(sum(cellData.loc[lab, 'CellGeneCount']), abs=0.05)
        assert 'soft' in cc['counts_are']
        got = {g['gene']: g['counts'] for g in cc['per_gene']}
        want = dict(zip(cellData.loc[lab, 'Genenames'], cellData.loc[lab, 'CellGeneCount']))
        for g, n in want.items():
            assert got[g] == pytest.approx(n, abs=0.002), (lab, g)


def test_cell_counts_for_one_gene(fitted):
    run, _, _ = fitted
    all_of_them = {g['gene']: g['counts'] for g in run.cell_counts(105)['per_gene']}
    for g, n in all_of_them.items():
        assert run.cell_counts(105, gene=g)['counts'] == pytest.approx(n)
    with pytest.raises(ValueError):
        run.cell_counts(105, gene='NoSuchGene')


def test_explain_cell_sums_are_the_sums_of_the_listed_diffs(fitted):
    """The two totals are the chart titles in the viewer's cell diagnostics. They
    are given so the agent does not add the diffs up itself, so they had better be
    exactly that sum."""
    run, _, _ = fitted
    for lab in LABELS:
        e = run.explain_cell(lab)
        got_a = sum(g['diff'] for g in e['genes_favouring_assigned'])
        got_o = sum(g['diff'] for g in e['genes_favouring_compared'])
        assert e['sum_favouring_assigned'] == pytest.approx(got_a, abs=1e-4), lab
        assert e['sum_favouring_compared'] == pytest.approx(got_o, abs=1e-4), lab
        assert e['sum_favouring_assigned'] > 0 >= e['sum_favouring_compared']


def test_explain_cell_counts_are_the_cell_counts(fitted):
    """The counts explain_cell quotes per gene must be the ones cell_counts gives
    for the same gene, so the two tools cannot show the user different numbers."""
    run, _, _ = fitted
    for lab in LABELS:
        e = run.explain_cell(lab)
        for g in e['genes_favouring_assigned'] + e['genes_favouring_compared']:
            assert g['counts'] == pytest.approx(run.cell_counts(lab, gene=g['gene'])['counts']), (lab, g['gene'])


# --------------------------------------------------------- spots_in_cell

def test_spots_in_cell_matches_geneData_and_says_it_is_hard(fitted):
    run, _, geneData = fitted
    for lab in LABELS:
        sc = run.spots_in_cell(lab)
        assert sc['total_spots'] == int((geneData.inside_cell == lab).sum())
        assert 'no probabilities' in sc['spots_are']
        assert all(isinstance(g['spots'], int) for g in sc['per_gene'])


def test_spots_in_cell_for_one_gene(fitted):
    run, _, geneData = fitted
    inside = geneData[geneData.inside_cell == 103]
    for g, n in inside.gene_name.value_counts().items():
        assert run.spots_in_cell(103, gene=g)['spots'] == n


def test_soft_and_hard_are_different_numbers(fitted):
    """If these ever agreed for every cell the two tools would be redundant."""
    run, _, _ = fitted
    soft = [run.cell_counts(l)['total_counts'] for l in LABELS]
    hard = [run.spots_in_cell(l)['total_spots'] for l in LABELS]
    assert any(abs(s - h) > 0.5 for s, h in zip(soft, hard))


# --------------------------------------------------------- the guard rails

def test_background_is_refused_with_a_clear_message(fitted):
    """The agent reads assigned_to: 'background' and may hand it straight back."""
    run, _, _ = fitted
    with pytest.raises(ValueError, match='background'):
        run.cell_counts('background')
    with pytest.raises(ValueError, match='background'):
        run.cell(0)


def test_unknown_cell_and_spot(fitted):
    run, _, _ = fitted
    with pytest.raises(KeyError):
        run.cell(999999)
    with pytest.raises(KeyError):
        run.explain_spot(999999999)


def test_open_run_on_an_empty_folder(tmp_path):
    with pytest.raises(FileNotFoundError):
        open_run(tmp_path)


# --------------------------------------------------------- spots_of_cell

def test_argmax_spots_match_geneData_neighbour(fitted):
    run, _, geneData = fitted
    for lab in LABELS:
        got = run.spots_of_cell(lab)
        assert got['definition'] == 'most likely parent'
        assert {r['spot'] for r in got['spots']} == set(geneData.index[geneData.neighbour == lab])
        assert 'not a hard assignment' in got['spots_are']


def test_thresholded_spots_match_geneData_and_sum_to_the_counts(fitted):
    run, cellData, geneData = fitted
    for lab in LABELS:
        got = run.spots_of_cell(lab, min_prob=0.0001)
        assert got['definition'] == 'prob > 0.0001'

        want = set()
        for sid, cands, probs in zip(geneData.index, geneData.neighbour_array, geneData.neighbour_prob):
            if lab in cands and probs[list(cands).index(lab)] > 0.0001:
                want.add(sid)
        assert {r['spot'] for r in got['spots']} == want
        assert got['sum_of_probs'] == pytest.approx(sum(cellData.loc[lab, 'CellGeneCount']), abs=0.05)
        assert 'decomposition of cell_counts' in got['spots_are']


def test_the_two_definitions_nest(fitted):
    """Argmax is a subset of any threshold below 0.5, and normally much smaller."""
    run, _, _ = fitted
    for lab in LABELS:
        argmax = {r['spot'] for r in run.spots_of_cell(lab)['spots']}
        wide = {r['spot'] for r in run.spots_of_cell(lab, min_prob=0.0001)['spots']}
        assert argmax <= wide


def test_spots_come_sorted_by_probability(fitted):
    run, _, _ = fitted
    probs = [r['prob'] for r in run.spots_of_cell(105, min_prob=0.0)['spots']]
    assert probs == sorted(probs, reverse=True)


def test_spots_of_cell_refuses_a_bad_label(fitted):
    run, _, _ = fitted
    with pytest.raises(KeyError):
        run.spots_of_cell(999999)


# ------------------------------------------ the tsv rows, value for value

@pytest.fixture(scope='module')
def tsv(fitted):
    """The two tsv files as the run wrote them, lists parsed back."""
    from pciSeq.src.core.io import read_tsv
    run, _, _ = fitted
    d = run.path.parent.parent.parent / 'tsv'
    return (read_tsv(str(d / 'cellData.tsv')).set_index('Cell_Num'),
            read_tsv(str(d / 'geneData.tsv')).set_index('spot_id'))


def test_cell_row_is_the_cellData_row(fitted, tsv):
    run, _, _ = fitted
    cells, _ = tsv
    for lab in LABELS:
        got = run.cell_row(lab)
        want = cells.loc[lab]
        assert got['Cell_Num'] == lab
        for col in ('X', 'Y', 'Z'):
            assert got[col] == pytest.approx(want[col], abs=1e-3), (lab, col)
        assert got['Genenames'] == list(want['Genenames'])
        assert got['CellGeneCount'] == pytest.approx(list(want['CellGeneCount']), abs=1e-3)
        assert got['ClassName'] == list(want['ClassName'])
        assert got['Prob'] == pytest.approx(list(want['Prob']), abs=1e-3)
        # spot_id: same spots under each gene, order inside a gene is not promised
        assert [sorted(s) for s in got['spot_id']] == [sorted(s) for s in want['spot_id']], lab


def test_cell_shape_is_the_covariance_the_run_fitted(fitted):
    """The shape tool reads the covariance diagnostics.db keeps; it must be the
    model's own cells.cov, and the radii its eigenvalues."""
    run, _, _ = fitted
    obj = pd.read_pickle(run.path.parent.parent.parent / 'debug' / 'pciSeq.pickle')
    for lab in LABELS:
        got = run.cell_shape(lab)
        want = obj.cells.cov[run.to_internal(lab)]
        assert np.allclose(got['cov'], want, atol=1e-3), lab
        radii = np.sqrt(np.sort(np.linalg.eigvalsh(want))[::-1])
        assert got['radii'] == pytest.approx(radii, abs=1e-3), lab
        assert got['radii'] == sorted(got['radii'], reverse=True)
        assert got['elongation'] == pytest.approx(radii[0] / radii[-1], abs=1e-3)
    with pytest.raises(ValueError):
        run.cell_shape(0)


def test_shape_words():
    from pciSeq.src.mcp.tools import shape_of
    # a sphere of the mean cell radius
    got = shape_of(100 * np.eye(3), mcr=10.0)
    assert got['shape'].startswith('round, the prior sphere')
    assert got['long_axis'] is None
    assert got['radii_vs_mean_cell_radius'] == [1.0, 1.0, 1.0]
    # twice as long along z as across
    got = shape_of(np.diag([25.0, 25.0, 100.0]), mcr=10.0)
    assert got['shape'] == 'elongated'
    assert got['radii'] == [10.0, 5.0, 5.0]
    assert got['elongation'] == 2.0
    assert got['long_axis'] == {'x': 0.0, 'y': 0.0, 'z': 1.0}
    assert 'mostly along z' in got['long_axis_is']
    # tilted 45 degrees in the xy plane, 2D
    got = shape_of(np.array([[50.0, 40.0], [40.0, 50.0]]))
    assert got['shape'] == 'elongated'
    assert 'between x and y' in got['long_axis_is']
    assert got['long_axis']['x'] == pytest.approx(0.707, abs=1e-3)
    assert 'radii_vs_mean_cell_radius' not in got


def test_spot_row_is_the_geneData_row(fitted, tsv):
    run, _, _ = fitted
    _, spots = tsv
    for sid in list(spots.index[:150]) + list(spots.index[-50:]):
        got = run.spot_row(sid)
        want = spots.loc[sid]
        assert got['gene_name'] == want['gene_name']
        assert got['gene_id'] == want['gene_id']
        for col in ('x', 'y', 'z', 'omp_score', 'omp_intensity'):
            assert got[col] == pytest.approx(want[col], abs=1e-3), (sid, col)
        assert got['plane_id'] == want['plane_id']
        assert got['neighbour'] == want['neighbour']
        assert got['neighbour_array'] == list(want['neighbour_array'])
        assert got['neighbour_prob'] == pytest.approx(list(want['neighbour_prob']), abs=1e-3)
        assert got['is_hard_misread'] == want['is_hard_misread']
        assert got['inside_cell'] == want['inside_cell']


def test_cell_row_spot_ids_are_the_soft_count_decomposition(fitted):
    """The probabilities of the listed spots add up to the counts, per gene."""
    run, _, _ = fitted
    row = run.cell_row(105)
    for gene, count, sids in zip(row['Genenames'], row['CellGeneCount'], row['spot_id']):
        total = sum(next(c['prob'] for c in run.explain_spot(s)['candidates'] if c['cell'] == 105)
                    for s in sids)
        assert total == pytest.approx(count, abs=0.002), gene


# --------------------------------------------------------- threads

def test_a_run_opened_on_one_thread_answers_on_another(fitted):
    """An MCP server runs each tool call on a worker thread and does not promise the
    same one twice, so open_run and the next call can land on different threads.
    sqlite refuses a connection used across threads unless told otherwise, and that
    showed up as a one-in-four failure of the full suite before it was caught.
    This forces the case every time."""
    import concurrent.futures as cf
    import threading

    run, cellData, _ = fitted
    here = threading.get_ident()
    with cf.ThreadPoolExecutor(max_workers=1) as pool:
        there = pool.submit(threading.get_ident).result()
        assert there != here
        got = pool.submit(run.cell_counts, 105).result()
        spot = pool.submit(run.explain_spot, 0).result()
    assert got['total_counts'] == pytest.approx(sum(cellData.loc[105, 'CellGeneCount']), abs=0.05)
    assert sum(c['prob'] for c in spot['candidates']) == pytest.approx(1.0, abs=1e-5)


def test_queries_from_many_threads_at_once(fitted):
    """The lock in Run._query. Hammer one connection from several threads and every
    answer has to be the right one, not a crash and not somebody else's row."""
    import concurrent.futures as cf

    run, _, _ = fitted
    want = {lab: run.cell(lab)['classes'][0]['class'] for lab in LABELS}
    with cf.ThreadPoolExecutor(max_workers=8) as pool:
        futures = [pool.submit(run.cell, lab) for lab in LABELS * 5]
        for lab, fut in zip(LABELS * 5, futures):
            assert fut.result()['classes'][0]['class'] == want[lab]


# --------------------------------------------------------- the narrative

def test_narrative_tells_the_story_without_units(fitted):
    """Every explain_cell answer carries a plain-words story built from its own
    numbers. It names both classes, the top gene and what settled it, and it never
    mentions a nat: the evidence comes out as odds or as a word, because the reader
    it is written for has never heard of a log-likelihood."""
    import re
    run, _, _ = fitted
    for lab in LABELS:
        e = run.explain_cell(lab)
        story = e['narrative']
        assert e['assigned'] in story and e['compared_with'] in story
        top = (e['genes_favouring_assigned'] or e['genes_favouring_compared'])[0]['gene']
        assert top in story
        assert re.search(r'So (the genes settled it|this call is the neighbourhood|'
                         r'the neighbourhood settled it|the prior settled it)', story), lab
        assert not re.search(r'\bnats?\b', story), lab
        # probabilities like 1.00 and 0.01 are allowed; a negative number or anything
        # 10 or more with a decimal point is a log-likelihood that leaked through
        assert not re.search(r'-\d|\b\d{2,}\.\d', story), 'raw score leaked: %d' % lab
        assert '\u2014' not in story and '\u2013' not in story


def test_narrative_when_the_neighbourhood_overrules_the_genes():
    """The fixture may not contain such a cell, so build the case by hand: genes
    favour the other class by 4 (about 55 to one), the neighbourhood favours the
    assigned one by 12. The story has to say the neighbourhood overruled the genes,
    name the class the genes alone would have picked, and give the odds in words."""
    from pciSeq.src.mcp.tools import narrate_cell
    e = {
        'cell': 4308, 'assigned': 'L5 ET', 'compared_with': 'L4/5 IT',
        'prob_assigned': 0.999, 'prob_compared': 0.001,
        'score': {'gene_loglik': {'assigned': -300.0, 'compared': -296.0},
                  'log_prior': {'assigned': -4.3, 'compared': -4.3},
                  'spatial': {'assigned': 12.5, 'compared': 0.5}},
        'genes_favouring_assigned': [{'gene': 'Cpne7', 'counts': 9.0, 'diff': 8.2}],
        'genes_favouring_compared': [{'gene': 'Car4', 'counts': 6.0, 'diff': -7.2}],
    }
    story = narrate_cell(e)
    assert 'looks more like L4/5 IT, about 55 to one, mostly because of Car4' in story
    assert 'The neighbourhood strongly favours L5 ET over L4/5 IT, which is what decides it' in story
    # the spatial term compares two classes, it says nothing of how many neighbours
    # are of either (espio cell 2413: 'mostly Ndnf' with 2 of 9 Ndnf)
    assert 'neighbouring cells are' not in story and 'mostly L5 ET' not in story
    assert 'neighbourhood overruling the genes' in story
    assert 'would have been called L4/5 IT' in story
    assert 'comes out for L5 ET, about 3,000 to one' in story


def test_narrative_when_the_genes_decide():
    from pciSeq.src.mcp.tools import narrate_cell
    e = {
        'cell': 1, 'assigned': 'A', 'compared_with': 'B',
        'prob_assigned': 1.0, 'prob_compared': 0.0,
        'score': {'gene_loglik': {'assigned': -100.0, 'compared': -133.6},
                  'log_prior': {'assigned': -4.3, 'compared': -4.3},
                  'spatial': {'assigned': 4.8, 'compared': 0.0}},
        'genes_favouring_assigned': [{'gene': 'Ndnf', 'counts': 8.0, 'diff': 15.3},
                                     {'gene': 'Rgs5', 'counts': 7.0, 'diff': 15.3}],
        'genes_favouring_compared': [{'gene': 'Npy', 'counts': 3.0, 'diff': -7.6}],
    }
    story = narrate_cell(e)
    assert 'The genes point to A, overwhelmingly, beyond any doubt.' in story
    assert 'The strongest evidence comes from Ndnf and Rgs5' in story
    assert 'at less than 0.01' in story
    assert 'So the genes settled it, and the neighbourhood agreed.' in story


def test_strength_wording():
    """The odds bands, so a change in rounding shows up here and not in a doc."""
    from pciSeq.src.mcp.tools import _strength
    assert _strength(0.2) == 'only just'
    assert _strength(0.8) == 'slightly'
    assert _strength(2.0) == 'about 7 to one'
    assert _strength(4.1) == 'about 60 to one'
    assert _strength(8.4) == 'about 4,400 to one'
    assert _strength(12.5) == 'overwhelmingly, beyond any doubt'


# ---------------------------------------------------- the spot narrative

def test_spot_narrative_tells_the_story_without_units(fitted):
    import re
    run, _, geneData = fitted
    for sid in list(geneData.index[:40]) + list(geneData.index[-10:]):
        e = run.explain_spot(sid)
        story = e['narrative']
        assert e['gene'] in story
        if e['assigned_to'] == 'background':
            assert 'assigned to the background' in story
        else:
            assert 'assigned to cell %d' % e['assigned_to'] in story
            assert re.search(r'decides it|the call is clear|a narrow one', story), sid
        assert not re.search(r'\bnats?\b', story)
        assert not re.search(r'-\d|\b\d{2,}\.\d', story), 'raw score leaked: %d' % sid
        assert '\u2014' not in story and '\u2013' not in story


def test_spot_candidates_carry_their_class(fitted):
    run, cellData, _ = fitted
    e = run.explain_spot(0)
    for c in e['candidates']:
        if c['cell'] != 'background':
            assert c['class'] == cellData.loc[c['cell'], 'ClassName'][0]


def _spot(assigned, cands, bg_prob, bg_misread=-14.0):
    rows = []
    for cell, cls, gauss, cx, sc, cg, prob in cands:
        rows.append({'cell': cell, 'class': cls, 'Gaussian fit': gauss, 'class expression': cx,
                     'cell scale': sc, 'cell-gene scale': cg, 'gene efficiency': -0.35,
                     'bonus': 0.0, 'sum': gauss + cx + sc + cg - 0.35, 'prob': prob})
    rows.append({'cell': 'background', 'misread': bg_misread, 'sum': bg_misread, 'prob': bg_prob})
    return {'spot': 7, 'gene': 'Synpr', 'position': {}, 'candidates': rows, 'assigned_to': assigned}


def test_spot_narrative_when_distance_carries_it():
    """The docs spot: nearest by a clear margin, every candidate the same class, the
    runner up holds more reads. Distance wins."""
    from pciSeq.src.mcp.tools import narrate_spot
    e = _spot(18223, [(18223, 'DG', -10.4, 0.94, -0.66, 0.18, 0.74),
                      (21574, 'DG', -12.8, 0.94, 0.10, 0.07, 0.13),
                      (17371, 'DG', -12.1, 0.94, -0.30, -0.54, 0.10)], 0.01)
    s = narrate_spot(e)
    assert 'assigned to cell 18223, fairly confidently, with probability 0.74' in s
    assert 'Cell 18223 is the nearest candidate: about 5 to one over cell 17371 and about 11 to one over cell 21574' in s
    assert 'Every candidate is a DG cell' in s
    assert 'cell 21574 holds more reads overall than its class predicts' in s
    assert 'distance decides it' in s


def test_spot_narrative_when_expression_carries_it():
    from pciSeq.src.mcp.tools import narrate_spot
    e = _spot(2, [(2, 'CA2', -11.0, 0.9, 0.5, 0.6, 0.66),
                  (1, 'CA1', -10.6, -1.2, 0.0, 0.0, 0.32)], 0.01)
    s = narrate_spot(e)
    assert 'Cell 2 is not the nearest: cell 1 is closer' in s
    assert 'Cell 2 is a CA2 cell, whose class expresses Synpr more than the class of cell 1, CA1' in s
    assert 'expression decides it' in s


def test_spot_narrative_when_the_background_wins():
    from pciSeq.src.mcp.tools import narrate_spot
    e = _spot('background', [(6482, 'VLMC', -15.6, -1.3, 0.6, -0.36, 0.02)], 0.95)
    s = narrate_spot(e)
    assert 'assigned to the background, with probability 0.95' in s
    assert 'takes it for a misread' in s
    assert 'cell 6482 is the nearest, and its class is VLMC' in s
    assert 'decides it' not in s


# ----------------------------------------- settings, convergence and planes

def test_run_info_reports_the_settings_and_how_it_ended(fitted):
    """The fixture is run with max_iter 40 and tolerance 0.001, so those two have to
    come straight back, along with a convergence verdict in words."""
    run, _, _ = fitted
    info = run.run_info()
    # the run date is the one thing everybody asks for first
    assert info['run_date'].endswith('Z') and info['run_date'].startswith('20')
    assert info['python_version'] and 'numpy' in info['package_versions']
    # the mean cell radius is written to the metadata now; it is what every cell's
    # Gaussian is initialised with, so it is positive and smaller than the image
    assert 0 < info['mean_cell_radius'] < 1000
    assert info['mean_cell_radius'] == float(run._meta('mcr'))
    assert info['settings']['max_iter'] == 40
    assert info['settings']['CellCallTolerance'] == 0.001
    assert info['voxel_size'] == [1, 1, 1] and info['is3D'] is False
    assert isinstance(info['converged'], bool) and info['iterations'] >= 1
    assert ('converged after' in info['ended']) == info['converged']
    assert 'label_map' not in info['settings']
    s = run.summary()
    assert s['has_settings'] and s['iterations'] == info['iterations']


def test_run_info_on_a_run_without_settings(fitted):
    """Runs written before the export carry no config row. The answer must say so
    rather than invent settings."""
    run, _, _ = fitted
    saved = run.config, run.run_record
    run.config, run.run_record = None, None
    try:
        info = run.run_info()
        assert info['settings'] is None and 'before pciSeq exported' in info['note']
        assert info['commit'] is not None
        pos = run.explain_spot(0)['position']
        assert 'plane' not in pos and 'no voxel_size' in pos['z_is']
    finally:
        run.config, run.run_record = saved


def test_explain_spot_gives_the_plane_on_an_anisotropic_3d_run(rng, tmp_path):
    """Three planes with voxel_size [1, 1, 4]: the model works in a z four times the
    plane index, and the tool has to hand back the plane, not the scaled z.

    The spots sit between planes, z_plane 1.6 and so on, as real spot callers give
    them. With whole planes only this test passed for weeks while diagnostics.db cut
    the z to a whole number and put one espio spot in ten on the plane below."""
    import pciSeq
    from scipy.sparse import coo_matrix
    from pciSeq.src.mcp.tools import open_run
    from tests.test_label_identifiers import _label_image, GENES, CLASSES

    spots = pd.DataFrame({
        'gene_name': rng.choice(GENES, 900),
        'x': rng.uniform(0, 79, 900).astype(np.float32),
        'y': rng.uniform(0, 79, 900).astype(np.float32),
        'z_plane': rng.uniform(0, 2, 900).astype(np.float32),
    })
    scref = pd.DataFrame(rng.random((len(GENES), len(CLASSES))) * 50, index=GENES, columns=CLASSES)
    scref.index.name = 'gene_name'
    lab = _label_image()
    # one matrix per plane: process_labels remaps each plane in place, so the same
    # object three times would be renumbered on the first pass and fail on the second
    cellData, geneData = pciSeq.fit(spots=spots, coo=[coo_matrix(lab) for _ in range(3)], scRNAseq=scref,
                                    opts={'max_iter': 3, 'save_data': True, 'output_path': str(tmp_path),
                                          'voxel_size': [1, 1, 4], 'CellCallTolerance': 0.5})
    run = open_run(tmp_path)
    assert run.run_info()['voxel_size'] == [1, 1, 4] and run.run_info()['is3D'] is True
    gd = geneData.set_index('spot_id')
    # the storage round trip: what the db gives back is the z the model had, and
    # the plane from it is the spot's plane_id
    assert (gd['z'] % 1 != 0).mean() > 0.9
    for sid in gd.index[:150]:
        pos = run.explain_spot(sid)['position']
        # geneData rounds z to 3 decimals, the db keeps it whole
        assert pos['z'] == pytest.approx(float(gd.loc[sid, 'z']), abs=5e-4), sid
        assert pos['plane'] == gd.loc[sid, 'plane_id'], sid
        assert run.spot(sid)['position']['plane'] == gd.loc[sid, 'plane_id'], sid

    # a file from before the fix: z cut to a whole number and no spot_z_exact. The
    # tools give no plane for it, and fix_spot_z puts the z back from arrow_spots
    import shutil
    import sqlite3
    from pciSeq.src.core.io import fix_spot_z
    db = run.path
    old_dir = tmp_path / 'old' / 'viewer_data' / 'diagnostics'
    old_dir.mkdir(parents=True)
    old_db = old_dir / 'diagnostics.db'
    shutil.copy(db, old_db)
    with sqlite3.connect(old_db) as con:
        con.execute('UPDATE spots SET z = CAST(z AS INTEGER)')
        con.execute("DELETE FROM metadata WHERE key = 'spot_z_exact'")
    spots_dir = db.parent.parent / 'arrow_spots'
    old_run = open_run(old_db)
    pos = old_run.explain_spot(int(gd.index[0]))['position']
    assert 'plane' not in pos and 'fix_spot_z' in pos['z_is']
    fix_spot_z(str(old_db), str(spots_dir))
    fixed = open_run(old_db)
    for sid in gd.index[:150]:
        pos = fixed.explain_spot(sid)['position']
        # arrow_spots keeps z to 3 decimals, as geneData does
        assert pos['z'] == pytest.approx(float(gd.loc[sid, 'z']), abs=1e-6), sid
        assert pos['plane'] == gd.loc[sid, 'plane_id'], sid
    fix_spot_z(str(old_db), str(spots_dir))     # a second time does nothing

    # the same cells on all three planes, so the segmentation puts every centroid at
    # plane 1 (z = 4 scaled). The fit then moves the centres a little in z, so the
    # expected counts come from the saved Z. find_cells(plane) has to compare planes,
    # not the scaled z, which is what floor(Z / 4) checks
    planes = np.floor(cellData['Z'].to_numpy() / 4).astype(int)
    assert (planes == 1).sum() > 0
    # cell_row gives the plane next to Z, so the agent never reads Z as a plane
    # (espio cell 8: Z 9.256 was answered as plane 9, it is plane 3)
    for lab, want in zip(cellData['Cell_Num'], planes):
        row = run.cell_row(int(lab))
        assert row['plane'] == want and 'not a plane' in row['Z_is'], lab
    assert run.find_cells(plane=1)['n_matching'] == (planes == 1).sum()
    assert run.find_cells(plane=0)['n_matching'] == (planes == 0).sum()
    assert run.find_cells(plane=4)['n_matching'] == 0


def test_plane_is_the_floor_like_the_spots_plane_id():
    """A centroid between two planes belongs to the lower one, as pciSeq gives a spot
    its plane_id and the viewer its Centroid Plane. Cell 2413 on espio is at plane
    21.56, the viewer says 21; the tools said 22 while they rounded."""
    from pciSeq.src.mcp.tools import Run
    r = Run.__new__(Run)
    r.config = {'voxel_size': [0.28, 0.28, 0.7]}
    assert r._plane_of(53.91) == 21
    # exactly on a plane, even when the ratio is not exact in floating point
    assert all(r._plane_of(p * 0.7 / 0.28) == p for p in range(90))


def test_narrative_says_when_a_gene_counts_by_absence():
    """Cell 16609 on espio: Pcp4 is the second strongest gene for CA3 with 0.00002
    counts, because a CA2 cell would hold some. Listed as 'the cell holds these in
    the amounts a CA3 cell does' that reads as a mistake, and the viewer's agent
    quietly dropped it. The story has to say absence is the evidence."""
    from pciSeq.src.mcp.tools import narrate_cell
    e = {
        'cell': 16609, 'assigned': 'CA3', 'compared_with': 'CA2',
        'prob_assigned': 1.0, 'prob_compared': 0.0,
        'score': {'gene_loglik': {'assigned': -200.0, 'compared': -220.0},
                  'log_prior': {'assigned': -4.3, 'compared': -4.3},
                  'spatial': {'assigned': 3.0, 'compared': 0.0}},
        'genes_favouring_assigned': [{'gene': 'Cdh9', 'counts': 6.1, 'diff': 9.8,
                                      'mean_in_assigned': 5.8, 'mean_in_compared': 0.3},
                                     {'gene': 'Pcp4', 'counts': 0.00002, 'diff': 2.55,
                                      'mean_in_assigned': 0.66, 'mean_in_compared': 3.6},
                                     {'gene': 'Tnfaip8l3', 'counts': 3.5, 'diff': 2.2,
                                      'mean_in_assigned': 3.1, 'mean_in_compared': 0.4}],
        'genes_favouring_compared': [{'gene': 'Amigo2', 'counts': 1.2, 'diff': -0.9,
                                      'mean_in_assigned': 0.5, 'mean_in_compared': 2.9},
                                     {'gene': 'Kcnq5', 'counts': 0.01, 'diff': -0.4,
                                      'mean_in_assigned': 1.25, 'mean_in_compared': 0.2}],
    }
    s = narrate_cell(e)
    assert 'The strongest evidence comes from Cdh9 and Tnfaip8l3' in s
    assert ('Pcp4 counts the same way by its absence: the cell holds almost none, where '
            'cells this run called CA2 hold about 3.6 on average.') in s
    assert 'A few genes, Amigo2, look more like CA2, but they are outweighed.' in s
    # 1.25 rounds half to even, 1.2, and the JS port must agree
    assert ('The near absence of Kcnq5 also leans towards CA2, since cells this run called '
            'CA3 hold about 1.2 on average, but not by enough.') in s

    # absence as the only evidence, two genes
    e['genes_favouring_assigned'] = [{'gene': 'Pcp4', 'counts': 0.0, 'diff': 2.55,
                                      'mean_in_assigned': 0.66, 'mean_in_compared': 3.6},
                                     {'gene': 'Car2', 'counts': 0.05, 'diff': 1.1,
                                      'mean_in_assigned': 0.1, 'mean_in_compared': 1.2}]
    s = narrate_cell(e)
    assert ('The strongest evidence is absence: the cell holds almost no Pcp4 and Car2, '
            'where cells this run called CA2 hold about 3.6 and 1.2 on average.') in s


def test_explain_cell_reports_the_class_means(fitted):
    """The Cell Inspector's Gene Expression table shows the average count over the
    cells this run called each class. explain_cell has to carry the same two numbers,
    or the story cannot back 'cells called CA2 hold about 3.6' with a figure."""
    run, _, _ = fitted
    means = np.asarray(run._meta('mean_gene_reads_per_class', parse=True), dtype=np.float32)
    names = list(run.class_names)
    e = run.explain_cell(101)
    rows = e['genes_favouring_assigned'] + e['genes_favouring_compared']
    assert rows
    for g in rows:
        gi = list(run.gene_panel).index(g['gene'])
        assert g['mean_in_assigned'] == pytest.approx(means[gi, names.index(e['assigned'])])
        assert g['mean_in_compared'] == pytest.approx(means[gi, names.index(e['compared_with'])])
    # not 'scRNAseq profile': the cell type definitions need not be single-cell data
    assert 'not the cell type definitions' in e['means_are']


# ------------------------------------------------------- the whole run

def test_gene_counts_add_up(fitted):
    """The soft counts of a gene split by class must sum to its counts in cells,
    and the hard split too; the top cells must really hold the most of it."""
    run, _, _ = fitted
    g = run.gene(str(run.gene_panel[0]))
    assert g['counts_in_cells'] == pytest.approx(sum(r['soft'] for r in g['counts_per_class']), abs=1e-2)
    assert g['counts_in_cells'] == pytest.approx(sum(r['in_cells_called'] for r in g['counts_per_class']), abs=1e-2)
    assert g['eta'] > 0 and g['eta_times_Inefficiency'] > 0 and g['total_spots'] > 0
    counts = [c['counts'] for c in g['top_cells']]
    assert counts == sorted(counts, reverse=True)
    assert g['top_cells'][0]['counts'] == pytest.approx(
        run.cell_counts(g['top_cells'][0]['cell'], gene=g['gene'])['counts'])
    with pytest.raises(ValueError):
        run.gene('NoSuchGene')


def test_theta_and_gamma_are_the_stored_arrays(fitted):
    run, _, _ = fitted
    c = run._cell(105)
    th = run.theta(105)
    assert th['theta'] == c['theta']
    best = th['theta_bar_per_class'][0]
    assert best['theta_bar'] == pytest.approx(c['theta_bar'][c['assigned_class_idx']])
    gm = run.gamma(105)
    assert len(gm['gamma']) == run.nG and gm['class'] == best['class']
    one = run.gamma(105, gene=gm['gamma'][3]['gene'])
    assert one['gamma'] == gm['gamma'][3]['gamma']
    assert 'assigned class only' in gm['gamma_is']
    # expected is the count the fit used: put back in the formula it gives gamma again
    r_spot = float(run._meta('rSpot'))
    for row in gm['gamma']:
        assert row['expected'] >= 0
        if row['expected'] > 0:
            assert (r_spot + row['counts']) / (r_spot + row['expected']) == pytest.approx(row['gamma'], rel=1e-4)
    assert one['expected'] == gm['gamma'][3]['expected']


def test_spot_agrees_with_explain_spot(fitted):
    run, _, geneData = fitted
    sid = geneData.index[5]
    s, e = run.spot(sid), run.explain_spot(sid)
    assert s['assigned_to'] == e['assigned_to'] and s['gene'] == e['gene']
    assert [c['prob'] for c in s['candidates']] == [c['prob'] for c in e['candidates']]
    with pytest.raises(KeyError):
        run.spot(-1)


def test_neighbours_match_the_cell_and_carry_distances(fitted):
    run, _, _ = fitted
    nb = run.neighbours(101)
    assert [r['cell'] for r in nb['neighbours']] == run.cell(101)['neighbours']
    assert all(r['distance_xy'] >= 0 for r in nb['neighbours'])
    assert nb['class'] == run.explain_cell(101)['assigned']


def test_class_counts_cover_every_cell(fitted):
    run, cellData, _ = fitted
    cc = run.class_counts()
    assert cc['n_cells'] == 16 == sum(r['cells'] for r in cc['classes'])
    assert sum(r['soft'] for r in cc['classes']) == pytest.approx(16, abs=1e-3)
    assert cc['classes'][0]['class'] == 'Zero'
    want = cellData['ClassName'].map(lambda v: v[0]).value_counts()
    for r in cc['classes']:
        assert r['cells'] == int(want.get(r['class'], 0)), r['class']
    assert run.class_counts(min_counts=1e9)['n_cells'] == 0


def test_class_prior_is_the_saved_log_prior(fitted):
    run, _, _ = fitted
    cp = run.class_prior()
    saved = dict(zip(run.class_names, run.metadata('log_prior')['value']))
    assert len(cp['classes']) == run.nK
    assert cp['classes'][0]['class'] == 'Zero'
    for r in cp['classes']:
        assert r['log_weight'] == pytest.approx(saved[r['class']])
        assert r['weight'] == pytest.approx(np.exp(saved[r['class']]))
    weights = [r['weight'] for r in cp['classes'][1:]]
    assert weights == sorted(weights, reverse=True)
    assert cp['cell_type_prior'] == run.config['cell_type_prior']


def test_find_cells_filters(fitted):
    run, _, _ = fitted
    e = run.explain_cell(101)
    hit = run.find_cells(class_name=e['assigned'])
    assert 101 in [r['cell'] for r in hit['cells']]
    assert hit['n_matching'] == run.class_counts()['classes'][
        [r['class'] for r in run.class_counts()['classes']].index(e['assigned'])]['cells']
    sure = run.find_cells(top_two_within=0.0)
    assert sure['n_matching'] <= 16
    assert run.find_cells(n=3)['shown'] == 3
    with pytest.raises(ValueError):
        run.find_cells(class_name='NoSuchClass')


def test_metadata_lists_and_returns(fitted):
    run, _, _ = fitted
    keys = {k['key']: k['kind'] for k in run.metadata()['keys']}
    assert keys['eta_bar'].startswith('list of %d' % run.nG)
    assert keys['mean_gene_reads_per_class'] == 'list of %d of lists of %d' % (run.nG, run.nK)
    assert run.metadata('nC')['value'] == run.nC
    assert run.metadata('class_names')['value'] == list(run.class_names)
    with pytest.raises(KeyError):
        run.metadata('no_such_key')


def test_spots_of_cell_by_gene(fitted):
    """The per gene count is the tool's, and it is the filter of the full list."""
    run, _, _ = fitted
    full = run.spots_of_cell(103)
    gene = full['spots'][0]['gene']
    mine = run.spots_of_cell(103, gene=gene)
    assert mine['n_spots'] == sum(1 for s in full['spots'] if s['gene'] == gene)
    assert mine['gene'] == gene
    with pytest.raises(ValueError):
        run.spots_of_cell(103, gene='NoSuchGene')


def test_the_odds_wording_survives_a_huge_difference():
    """Cell 2413 on espio against 327 Oligo NN: the gene difference is far past
    e**709 and math.exp raised OverflowError, so explain_cell crashed on a plain
    'why is it not an oligodendrocyte'."""
    from pciSeq.src.mcp.tools import _strength
    assert _strength(5000.0) == 'overwhelmingly, beyond any doubt'
    assert _strength(1e9) == 'overwhelmingly, beyond any doubt'
    assert _strength(3.0) == 'about 20 to one'


def test_spot_narrative_compares_the_two_classes():
    """Spot 18377 on espio: both classes express Cnr1, both class expression terms
    are below zero, and the narrative said one class 'barely expresses' it and the
    other 'does not', the wrong way round. The term is a log, its sign means nothing
    on its own, so the narrative compares the two cells."""
    from pciSeq.src.mcp.tools import narrate_spot
    e = _spot(12693, [(12693, '005 L5 IT CTX Glut', -12.65, -0.62, 0.39, 0.31, 0.48),
                      (9077, '022 L5 ET CTX Glut', -13.84, -0.34, 0.71, 0.86, 0.47)], 0.04)
    s = narrate_spot(e)
    assert 'whose class expresses Synpr less than the class of cell 9077, 022 L5 ET CTX Glut' in s
    assert 'does not' not in s and 'barely' not in s
    e = _spot(12693, [(12693, 'A', -12.65, -0.62, 0.39, 0.31, 0.48),
                      (9077, 'B', -13.84, -0.55, 0.71, 0.86, 0.47)], 0.04)
    assert 'the two classes express Synpr about equally' in narrate_spot(e)

