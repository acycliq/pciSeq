"""Answer questions about a finished run.

These are the tools an agent calls. They are plain functions on purpose: the logic
should not know whether it is being reached over MCP, from a viewer, or from a
notebook. The MCP server is a thin wrapper on top.

Everything comes out of diagnostics.db, which is self contained: its metadata table
carries the gene panel, the class names, eta, the class means, the prior and the
label map, so no pickle and none of the big tsv files are needed.

Three rules hold for every tool here:

1. cell labels in and out are the labels of your segmentation, never the internal
   ones. See `to_internal` in cell_utils for why the two differ.
2. counts are soft, weighted by the assignment probability, and every answer that
   reports one says so. Containment is hard and says that too.
3. row 0 of the model arrays is the background pseudocell and is never reported as
   a cell.

One thing to know about precision. The saved files keep 3 decimals, so a probability
under 0.0005 reads as 0.0 in them. cell_row and spot_row hand back exactly what the
file says, which is right for "what does the file report", but it means neither can
show you a probability finer than that. explain_spot recomputes from diagnostics.db
and gives the full float32, so use that when the small numbers matter. Note also that
cellData truncates and geneData rounds, summary.py:31 against summary.py:111, so the
last digit of the two files is not arrived at the same way.
"""
import ast
import json
import logging
import math
import operator
import sqlite3
import threading
from pathlib import Path

import numpy as np

from ..core.utils.geometry import plane_of

logger = logging.getLogger(__name__)

_MISSING = object()

# counts below this are dropped from the per gene listings. Same threshold cellData
# uses, so the two agree about which genes a cell holds.
_COUNT_TOL = 0.001

# what a cell's blob columns unpack to, and the shape to give them. scaled_means,
# area_factor * mean_expression_adj before eta and theta, is a blob in format 1
# and rebuilt from the cell's area_factor and the shared matrix in format 2.
_CELL_BLOBS = {
    'theta_bar': ('nK',),
    'gene_count': ('nG',),
    'class_prob': ('nK',),
    'gamma_assigned': ('nG',),      # the winning class only, see the class docstring
    'mrf': ('nK',),
}


class _CellIndex:
    """The spot files' candidate lists turned inside out: for each cell, its
    entries. Built by Run._cell_index from the arrow spot shards."""

    def __init__(self, shards):
        import pyarrow.compute as pc
        spot_id, gene_id, cells, probs, first, parent = [], [], [], [], [], []
        base = 0
        gene_of = {}
        for tb, gmap in shards:
            gene_of = gmap
            n = tb.num_rows
            spot_id.append(tb['spot_id'].to_numpy())
            gene_id.append(tb['gene_id'].to_numpy())
            na = tb['neighbour_array'].combine_chunks()
            npb = tb['neighbour_prob'].combine_chunks()
            if n == 0 or len(na.values) == 0:
                base += n
                continue
            par = pc.list_parent_indices(na).to_numpy(zero_copy_only=False)
            offs = na.offsets.to_numpy()
            at = np.arange(len(par), dtype=np.int64) + offs[0]
            cells.append(pc.list_flatten(na).to_numpy(zero_copy_only=False).astype(np.int64))
            probs.append(pc.list_flatten(npb).to_numpy(zero_copy_only=False).astype(np.float32))
            first.append(at == offs[par])
            parent.append((par + base).astype(np.int32))
            base += n
        self.gene_of = gene_of
        self.spot_id = np.concatenate(spot_id) if spot_id else np.zeros(0, np.int64)
        self.gene_id = np.concatenate(gene_id) if gene_id else np.zeros(0, np.int64)
        cells = np.concatenate(cells) if cells else np.zeros(0, np.int64)
        # stable, so within a cell the entries keep the order of the files
        order = np.argsort(cells, kind='stable')
        cells = cells[order]
        self.labels, self.start = np.unique(cells, return_index=True)
        self.stop = np.append(self.start[1:], len(cells))
        self.spot_row = (np.concatenate(parent) if parent else np.zeros(0, np.int32))[order]
        self.prob = (np.concatenate(probs) if probs else np.zeros(0, np.float32))[order]
        self.first = (np.concatenate(first) if first else np.zeros(0, bool))[order]

    def entries_of(self, label):
        """(spot rows, probabilities, is-first flags) of the entries naming this
        cell, empty when no spot does."""
        i = np.searchsorted(self.labels, label)
        if i == len(self.labels) or self.labels[i] != label:
            return (np.zeros(0, np.int32), np.zeros(0, np.float32), np.zeros(0, bool))
        a, b = self.start[i], self.stop[i]
        return self.spot_row[a:b], self.prob[a:b], self.first[a:b]


class Run:
    """A finished pciSeq run, opened for questions.

    Parameters
    ----------
    path : str or Path
        The run folder, or the diagnostics.db itself. Anything above the database
        works, the file is searched for underneath.

    Notes
    -----
    diagnostics.db keeps `gamma_assigned`, the gamma of the winning class only, not
    the whole (nC, nG, nK) array. Questions about gamma under a different class need
    the pickle.
    """

    def __init__(self, path):
        self.path = _find_db(path)
        # check_same_thread=False because an MCP server runs each tool call on a
        # worker thread, and the pool does not promise the same thread twice: open_run
        # would build this connection on one thread and cell_counts hit it from
        # another, which sqlite refuses by default. It only started showing up as a
        # one-in-four failure in the full test suite. sqlite3.threadsafety is 1 in
        # the usual build, meaning the module is safe but a shared connection is not,
        # so every query also goes through _query, which holds a lock.
        self._con = sqlite3.connect('file:%s?mode=ro' % self.path, uri=True,
                                    check_same_thread=False)
        self._lock = threading.Lock()
        self._meta_cache = {}
        self._docs_pages = _MISSING
        self._mea = None            # mean_expression_adj, read once (format 2)

        self.nC = int(self._meta('nC'))
        self.nG = int(self._meta('nG'))
        self.nK = int(self._meta('nK'))
        self.nS = int(self._meta('nS'))
        self.gene_panel = np.array(self._meta('gene_panel', parse=True))
        self.class_names = np.array(self._meta('class_names', parse=True))

        # label_map comes back from JSON with string keys
        raw = self._meta('label_map', parse=True, default=None)
        self.label_map = {int(k): int(v) for k, v in raw.items()} if raw else None
        self._reverse_map = ({v: k for k, v in self.label_map.items()}
                             if self.label_map else None)

        # the key column of the cells table. It was cell_id before September 2026,
        # which clashed with cell_id in the other outputs, where it is the
        # segmentation label. Here it has always been the internal label.
        # old runs: remove the cell_id fallback with the next version bump, see cz1.3
        with self._lock:
            cols = [r[1] for r in self._con.execute('PRAGMA table_info(cells)')]
        self._cell_key = 'internal_label' if 'internal_label' in cols else 'cell_id'
        self._has_nbrs = 'neighbours' in cols
        # the covariance of each cell's Gaussian, saved since October 2026
        self._has_cov = 'cov' in cols

        # the resolved settings and the convergence record. Runs written before
        # September 2026 have neither; the tools say so rather than guess.
        # The file's format version. 0 means a run from before the key existed
        # (pre-1, to be refused once the old-run grace period ends, see cz1.9.1);
        # bigger than 1 means a future pciSeq broke the only-add promise, and
        # guessing at such a file risks quietly wrong answers, so refuse.
        self.format_version = int(self._meta('format_version', default='0'))
        if self.format_version > 2:
            raise NotImplementedError(
                'this run has diagnostics.db format version %d, and this pciSeq only '
                'knows versions 1 and 2. Read it with the pciSeq that made it, or newer.'
                % self.format_version)
        # whether the spots z has its decimals. Files from before October 2026 cut
        # it to a whole number, and a plane worked out from that is one too low for
        # about one spot in ten, so those files give no plane for a spot
        self.spot_z_exact = self._meta('spot_z_exact', default=None) == '1'
        # format 2 (October 2026): area_factor per cell in place of the scaled_means
        # blob, the spot blobs without their background slot, the candidate ids as
        # int32, no gene_inefficiency column (logeta_bar in the metadata instead) and
        # no bonus column on a run without InsideCellBonus. _cell and _spot_row hand
        # back the same arrays for either version.
        self._has_area_factor = 'area_factor' in cols
        # the run's own class score (gene_loglik per class, the per gene numbers of
        # the assigned class and the runner up), saved since 28 Sep 2026, and the
        # counts that score was computed from, since 15 Oct 2026
        self._has_saved_score = 'gene_loglik' in cols
        self._has_scored_counts = 'scored_gene_count' in cols
        with self._lock:
            self._spot_cols = {r[1] for r in self._con.execute('PRAGMA table_info(spots)')}

        self.config = self._meta('config', parse=True, default=None)
        self.run_record = self._meta('run', parse=True, default=None)

    # ------------------------------------------------------------ plumbing

    def _query(self, sql, params=()):
        """One row from the db, under the lock. See __init__ for why the lock."""
        with self._lock:
            return self._con.execute(sql, params).fetchone()

    def docs_pages(self):
        """The documentation pages saved inside this run (page -> markdown), or None
        for a run made before they were saved. Same table the viewer reads."""
        if self._docs_pages is _MISSING:
            with self._lock:
                has = self._con.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name='docs'").fetchone()
                rows = self._con.execute('SELECT page, text FROM docs').fetchall() if has else []
            self._docs_pages = {r[0]: r[1] for r in rows} or None
        return self._docs_pages

    def _meta(self, key, parse=False, default=_MISSING):
        """One row of the metadata table.

        `parse` decodes the value as JSON. `default` is returned when the key is
        absent, otherwise a missing key raises.
        """
        if key not in self._meta_cache:
            row = self._query(
                "SELECT value FROM metadata WHERE key = ?", (key,))
            self._meta_cache[key] = row[0] if row else None
        val = self._meta_cache[key]
        if val is None:
            if default is _MISSING:
                raise KeyError('diagnostics.db has no metadata key %r' % key)
            return default
        return json.loads(val) if parse else val

    def to_internal(self, label):
        """Segmentation label to the row it occupies in the arrays."""
        try:
            label = int(label)
        except (TypeError, ValueError):
            raise ValueError('%r is not a cell label. Spots assigned to the background '
                             'belong to no cell, so there is nothing to look up.' % (label,))
        if self.label_map is None:
            return label
        try:
            return self.label_map[label]
        except KeyError:
            raise KeyError('no cell %s in this segmentation' % label)

    def to_external(self, row):
        """The row of an array back to the segmentation label."""
        if self._reverse_map is None:
            return int(row)
        try:
            return self._reverse_map[int(row)]
        except KeyError:
            raise KeyError('no cell at row %s, this run has %d' % (row, self.nC))

    def _cell(self, label):
        """Every array of one cell, keyed by column name, taking a segmentation label."""
        row = self.to_internal(label)
        if row == 0:
            raise ValueError('cell 0 is the background pseudocell, not a cell')
        scaled = 'area_factor' if self._has_area_factor else 'scaled_means'
        cols = list(_CELL_BLOBS) + ['theta', 'assigned_class_idx', scaled]
        got = self._query(
            "SELECT %s FROM cells WHERE %s = ?" % (', '.join(cols), self._cell_key), (row,))
        if got is None:
            raise KeyError('cell %s is not in diagnostics.db' % label)

        out = {}
        sizes = {'nG': self.nG, 'nK': self.nK}
        for name, shape in _CELL_BLOBS.items():
            arr = np.frombuffer(got[cols.index(name)], dtype=np.float32)
            out[name] = arr.reshape([sizes[s] for s in shape])
        if self._has_area_factor:
            out['scaled_means'] = (np.float32(got[cols.index(scaled)])
                                   * self._mean_expression_adj()).astype(np.float32)
        else:
            out['scaled_means'] = np.frombuffer(got[cols.index(scaled)],
                                                dtype=np.float32).reshape(self.nG, self.nK)
        out['theta'] = float(got[cols.index('theta')])
        out['assigned_class_idx'] = int(got[cols.index('assigned_class_idx')])
        out['row'] = row
        if self._has_scored_counts:
            b = self._query('SELECT scored_gene_count FROM cells WHERE %s = ?' % self._cell_key, (row,))[0]
            out['scored_gene_count'] = np.frombuffer(b, dtype=np.float32)
        return out

    def _saved_score(self, row):
        """The class score the run saved for one cell, or None on a run from before
        it did: gene_loglik per class, runner_up_idx, and the per gene numbers of the
        assigned class and the runner up."""
        if not self._has_saved_score:
            return None
        gl, ru, ca, cr = self._query('SELECT gene_loglik, runner_up_idx, contr_assigned, '
                                     'contr_runner_up FROM cells WHERE %s = ?' % self._cell_key, (row,))
        f = lambda b: np.frombuffer(b, dtype=np.float32).astype(np.float64)
        return {'gene_loglik': f(gl), 'runner_up_idx': int(ru),
                'contr_assigned': f(ca), 'contr_runner_up': f(cr)}

    def _checked_contributions(self, c, saved, assigned, other):
        """The per gene numbers of the assigned class and another class, rebuilt,
        or None when the rebuild is not the run's own.

        The run saves per gene numbers for the runner up only, so for any other class
        they are rebuilt (_nb_contributions). They are used only when the rebuild of
        the assigned class matches the saved per gene numbers gene by gene, and the
        rebuilt totals of both classes match the saved gene_loglik: then the formula
        and the inputs are the run's. A change to the model the rebuild does not
        follow, or a run from before the scored counts were saved where the counts
        moved, fails here, and the caller says the detail is not available rather than
        show numbers the run would not have given. The viewer's nbContributions
        applies the same rule with the same tolerances."""
        contr = self._nb_contributions(c).astype(np.float64)
        ok = lambda r, s, rel, abs_: bool(np.all(np.abs(r - s) <= rel * np.abs(s) + abs_))
        if not ok(contr[:, assigned], saved['contr_assigned'], 1e-4, 1e-3):
            return None
        for k in (assigned, other):
            if not ok(contr[:, k].sum(), saved['gene_loglik'][k], 1e-4, 0.05):
                return None
        return contr[:, assigned], contr[:, other]

    # --------------------------------------------------------------- tools

    def commit(self):
        """The commit of the pciSeq that made this run, or None for a run from before
        it was recorded (the source tools refuse those rather than guess)."""
        return self._meta('pciSeq_provenance', parse=True, default={}).get('commit') or None

    def summary(self):
        """What this run is: its size, the code that made it, and what it can answer."""
        prov = self._meta('pciSeq_provenance', parse=True, default={})
        return {
            'path': str(self.path),
            'format_version': self.format_version,
            'cells': self.nC - 1,          # row 0 is the background
            'spots': self.nS,
            'genes': self.nG,
            'classes': self.nK,
            'pciSeq_version': prov.get('version'),
            'commit': prov.get('commit'),
            'labels_were_renumbered': self.label_map is not None,
            'has_containment': self._has_containment(),
            'has_settings': self.config is not None,
            'iterations': self.run_record['iterations'] if self.run_record else None,
            'converged': self.run_record['converged'] if self.run_record else None,
        }

    def run_info(self):
        """What produced this run and how it ended: the pciSeq version and commit, the
        resolved settings, the number of iterations and whether the loop converged.

        The provenance (version, commit, its date, when the run was made, python
        and the package versions) is always there. Runs written before the
        settings were exported carry nothing else; the answer says so.
        """
        prov = self._meta('pciSeq_provenance', parse=True, default={})
        out = {
            'path': str(self.path),
            'format_version': self.format_version,
            'pciSeq_version': prov.get('version'),
            'commit': prov.get('commit'),
            'commit_date': prov.get('commit_date'),
            'branch': prov.get('branch'),
            'run_date': prov.get('created_at'),
            'python_version': prov.get('python_version'),
            'os': prov.get('os'),
            'package_versions': prov.get('package_versions'),
            'cells': self.nC - 1, 'spots': self.nS, 'genes': self.nG, 'classes': self.nK,
            'mean_cell_radius': self.mean_cell_radius(),
            'mean_cell_radius_is': 'in pixels of the xy plane, the mean over the segmented '
                                   'cells of sqrt(area / pi), halved; unless cell_radius is set, '
                                   'the radius of the prior on every cell\'s Gaussian, whose '
                                   'covariance is then fitted per cell',
        }
        if self.config is None:
            out['settings'] = None
            out['note'] = ('this run was written before pciSeq exported its settings and '
                           'convergence record to diagnostics.db, so only the provenance '
                           'above is known. Rerunning with the current pciSeq records them.')
            return out
        cfg = self.config
        out['settings'] = cfg
        out['is3D'] = cfg.get('is3D')
        out['voxel_size'] = cfg.get('voxel_size')
        # iterations is None when the run was exported without the loop ever
        # finishing an iteration; there is no convergence story to tell then
        if self.run_record and self.run_record.get('iterations') is not None:
            out['iterations'] = self.run_record['iterations']
            out['converged'] = self.run_record['converged']
            delta = self.run_record.get('delta') or []
            out['final_delta'] = delta[-1] if delta else None
            out['tolerance'] = cfg.get('CellCallTolerance')
            out['ended'] = ('converged after %d iterations, the largest change in a spot '
                            'assignment fell below %s' % (out['iterations'], out['tolerance'])
                            if out['converged'] else
                            'stopped at max_iter, %d iterations, without converging: the '
                            'largest change was still %.4f against a tolerance of %s'
                            % (out['iterations'], out['final_delta'] or float('nan'),
                               out['tolerance']))
        return out

    def mean_cell_radius(self):
        """The mean cell radius, mcr, in pixels of the xy plane.

        Written to the metadata since September 2026. Older runs have it only in
        disguise: the Gaussian update is not run, so every cell keeps its initial
        covariance mcr^2 I, and cellData's sphere_scale is 3 sqrt of its eigenvalues,
        so sphere_scale / 3 on any row is mcr. None when neither is there.
        """
        val = self._meta('mcr', default=None)
        if val is not None:
            return float(val)
        row = self._tsv_row('cellData.tsv', r'^\d+\t')
        if row and isinstance(row.get('sphere_scale'), (list, tuple)) and row['sphere_scale']:
            return float(row['sphere_scale'][0]) / 3.0
        return None

    def _plane_of(self, z):
        """The plane index a scaled z belongs to, or None without voxel_size.

        Floor, not round, see plane_of in core/utils/geometry.py for why.
        """
        if not self.config or not self.config.get('voxel_size'):
            return None
        return plane_of(z, self.config['voxel_size'])

    def _spot_position(self, got):
        """x, y, z and the plane of a spot row, the plane only when it can be trusted."""
        plane, z_is = self._spot_plane(got[3])
        out = {'x': got[1], 'y': got[2], 'z': got[3]}
        if plane is not None:
            out['plane'] = plane
        out['z_is'] = z_is
        return out

    def _spot_plane(self, z):
        """The plane of a spot from its stored z, or None and the reason why not."""
        if not self.spot_z_exact:
            return None, ('the anisotropy scaled z the model works in, cut to a whole '
                          'number by the pciSeq that saved this run, so the plane cannot '
                          'be given. pciSeq.src.core.io.fix_spot_z repairs the file')
        plane = self._plane_of(z)
        if plane is None:
            return None, ('the anisotropy scaled z the model works in, not the plane '
                          'index. This run carries no voxel_size, so the plane cannot '
                          'be given')
        return plane, 'the anisotropy scaled z the model works in; plane is the plane index it sits on'

    def cell(self, label):
        """The headline facts about one cell."""
        c = self._cell(label)
        # stable: probabilities tie at exactly 0 and 1 in float32, and an unstable
        # sort would list the tied ones in an arbitrary order
        order = np.argsort(-c['class_prob'], kind='stable')
        counts = c['gene_count']
        top_genes = np.argsort(-counts, kind='stable')[:10]
        return {
            'cell': int(label),
            'total_counts': float(counts.sum()),
            'counts_are': 'soft, weighted by the spot assignment probabilities',
            'classes': [{'class': str(self.class_names[k]), 'prob': float(c['class_prob'][k])}
                        for k in order[:5]],
            'top_genes': [{'gene': str(self.gene_panel[g]), 'counts': float(counts[g])}
                          for g in top_genes if counts[g] > _COUNT_TOL],
            'theta': c['theta'],
            'theta_is': 'the soft scalar, averaged over the class probabilities, '
                        'not theta_bar of the assigned class',
            # None on runs written before diagnostics.db kept them
            'neighbours': self._neighbours(label),
            'neighbours_are': 'the cells the spatial (mrf) term listens to, nearest first',
        }

    def cell_shape(self, label):
        """The fitted shape of one cell: round or elongated, by how much, and which way.

        The covariance of the cell's Gaussian, read from diagnostics.db, as the
        ellipsoid it describes: the standard deviation along each of its own axes,
        longest first, the ratio of the longest to the shortest, and the direction of
        the longest. Compared with the mean cell radius, the radius of the prior
        sphere. Runs saved before October 2026 do not keep the covariance and the
        answer says so.
        """
        label = int(label)
        row = self.to_internal(label)
        if row == 0:
            raise ValueError('cell 0 is the background pseudocell, not a cell')
        if not self._has_cov:
            return {'cell': label, 'available': False,
                    'why': 'this run was saved before diagnostics.db kept the covariance '
                           'of each cell (October 2026), so its fitted shape is not '
                           'available'}
        got = self._query("SELECT cov FROM cells WHERE %s = ?" % self._cell_key, (row,))
        if got is None:
            raise KeyError('cell %s is not in diagnostics.db' % label)
        flat = np.frombuffer(got[0], dtype=np.float32)
        d = int(round(np.sqrt(flat.size)))
        out = {'cell': label}
        out.update(shape_of(flat.reshape(d, d), self.mean_cell_radius()))
        return out

    def explain_cell(self, label, vs_class=None, top_n=10):
        """Why this cell got its class, gene by gene.

        The class score is the sum of the per gene negative binomial contributions,
        plus the class prior, plus the spatial term. This recomputes those from the
        stored arrays and compares the assigned class against another one.

        The recomputation uses the final eta, and the loop updates eta after the last
        class update, so it reproduces the stored probabilities exactly only once the
        run has converged. On an unconverged run the argmax still agrees but the
        probabilities can be off in the second decimal.

        Parameters
        ----------
        label : int
            Cell label, as in your segmentation.
        vs_class : str, optional
            The class to compare against. Defaults to the runner up.
        top_n : int, default 10
            How many genes to report on each side.
        """
        c = self._cell(label)
        assigned = c['assigned_class_idx']
        if vs_class is None:
            order = np.argsort(-c['class_prob'])
            other = int(order[1] if order[0] == assigned else order[0])
            # the runner up is the second best class however small, but when no other
            # class has any probability at all there is none: argsort would just hand
            # back the first class of the list and the story would call it 'the
            # closest alternative'. Ask instead, the same rule as the viewer's chat.
            if len(order) < 2 or not c['class_prob'][other] > 0:
                raise ValueError('cell %s is %s and no other class has any probability, so '
                                 'there is no runner up to compare against. Ask the user '
                                 'which class to compare against and call again with '
                                 'vs_class' % (label, self.class_names[assigned]))
        else:
            hit = np.where(self.class_names == vs_class)[0]
            if not len(hit):
                raise ValueError('no class %r in this run' % vs_class)
            other = int(hit[0])
            if other == assigned:
                raise ValueError('cell %s is already %s, pick another class to compare'
                                 % (label, vs_class))

        # the per gene numbers: the run's own against the runner up, a checked rebuild
        # against any other class, and on a run from before the score was saved the
        # rebuild alone, as it always was
        saved = self._saved_score(c['row'])
        detail_note = None
        if saved is None:
            contr = self._nb_contributions(c)
            ca, co = contr[:, assigned], contr[:, other]
            gl_a, gl_o = float(ca.sum()), float(co.sum())
        else:
            gl_a, gl_o = float(saved['gene_loglik'][assigned]), float(saved['gene_loglik'][other])
            if other == saved['runner_up_idx']:
                ca, co = saved['contr_assigned'], saved['contr_runner_up']
            else:
                got = self._checked_contributions(c, saved, assigned, other)
                if got is None:
                    ca = co = None
                    detail_note = ('the per gene contributions against %s could not be '
                                   'rebuilt to match the run\'s own gene score, so only the '
                                   'three score totals are available. The gene log-likelihood '
                                   'difference above is still the run\'s own number.'
                                   % self.class_names[other])
                else:
                    ca, co = got
        diff = (ca - co) if ca is not None else np.zeros(self.nG)
        favours_assigned = np.argsort(-diff)[:top_n]
        favours_other = np.argsort(diff)[:top_n]
        log_prior = np.asarray(self._meta('log_prior', parse=True), dtype=np.float32)
        # (nG, nK): the average count of each gene over this run's cells, weighted by
        # their probability of being each class. The same numbers as the Gene
        # Expression table of the viewer's Cell Inspector, and the figure behind 'a
        # CA2 cell would hold some'. Not what the likelihood scores against, that is
        # the scRNAseq profile times the scale factors; this is a summary of the run.
        means = np.asarray(self._meta('mean_gene_reads_per_class', parse=True), dtype=np.float32)

        def side(idx, keep):
            return [{'gene': str(self.gene_panel[g]),
                     'counts': float(c['gene_count'][g]),
                     'mean_in_assigned': float(means[g, assigned]),
                     'mean_in_compared': float(means[g, other]),
                     'diff': float(diff[g])}
                    for g in idx if keep(diff[g])]

        for_a = side(favours_assigned, lambda d: d > 0) if ca is not None else []
        for_o = side(favours_other, lambda d: d < 0) if ca is not None else []
        total_a = gl_a + float(log_prior[assigned]) + float(c['mrf'][assigned])
        total_o = gl_o + float(log_prior[other]) + float(c['mrf'][other])
        out = {
            'cell': int(label),
            'assigned': str(self.class_names[assigned]),
            'compared_with': str(self.class_names[other]),
            'prob_assigned': float(c['class_prob'][assigned]),
            'prob_compared': float(c['class_prob'][other]),
            'score': {
                'gene_loglik': {'assigned': gl_a, 'compared': gl_o},
                'log_prior': {'assigned': float(log_prior[assigned]),
                              'compared': float(log_prior[other])},
                'spatial': {'assigned': float(c['mrf'][assigned]),
                            'compared': float(c['mrf'][other])},
                # the class score is the sum of the three parts; added here so the
                # agent never has to, the same as the viewer's explain_cell
                'total': {'assigned': total_a, 'compared': total_o},
                'difference': total_a - total_o,
                'difference_is': 'total assigned minus total compared, the three parts '
                                 'already added up; it equals log(prob_assigned / '
                                 'prob_compared). Quote it, never add the parts yourself',
            },
            'genes_favouring_assigned': for_a,
            'genes_favouring_compared': for_o,
            # the totals of the two lists, the sums in the chart titles of the viewer's
            # cell diagnostics. Given so the agent quotes them rather than adds up the
            # diffs itself, which it does badly.
            'sum_favouring_assigned': float(sum(g['diff'] for g in for_a)),
            'sum_favouring_compared': float(sum(g['diff'] for g in for_o)),
            'shared_genes': self._shared_genes(c, diff, means, assigned, other) if ca is not None else [],
            'shared_genes_are': 'the genes the cell holds most of that barely separate the two '
                                'classes, each fitting both within a factor of e. Both classes '
                                'express them, which is why these two were the finalists; the '
                                'call rests on the genes in the two lists above',
            'counts_are': 'soft, weighted by the spot assignment probabilities',
            'means_are': 'the average count over the cells of this run, each weighted by '
                         'its probability of being that class. An after the fact summary '
                         'of the run, not the cell type definitions the likelihood scores '
                         'against',
        }
        if detail_note:
            # no per gene detail, so no story to tell from it: the totals and the
            # note, the same as the viewer's explain_cell
            out['detail_note'] = detail_note
            return out
        out['narrative'] = narrate_cell(out)
        return out

    def _shared_genes(self, c, diff, means, assigned, other, top=10, n=5, within=1.0):
        """The cell's biggest counts that barely tell the two classes apart.

        Of the `top` genes the cell holds most of, the ones whose log-likelihood
        difference between the two classes is under `within` in size, biggest count
        first, at most `n`. These are the genes a model would otherwise mention with
        'both classes express Cnr1 and Htr3a', from its own memory; this gives it the
        data instead.
        """
        counts = c['gene_count']
        biggest = [g for g in np.argsort(-counts)[:top] if counts[g] > _COUNT_TOL]
        return [{'gene': str(self.gene_panel[g]),
                 'counts': float(counts[g]),
                 'mean_in_assigned': float(means[g, assigned]),
                 'mean_in_compared': float(means[g, other]),
                 'diff': float(diff[g])}
                for g in biggest if abs(diff[g]) < within][:n]

    def explain_spot(self, spot_id):
        """Why this spot went to the cell it did, term by term.

        One row per candidate cell plus the background. The terms are the ones the
        model used in its last spot update, so the probabilities here are the ones
        in geneData, at full float32 precision rather than the 3 decimals the file
        keeps. Use this rather than spot_row when a probability below 0.0005 matters.
        """
        # sqlite3 quietly matches nothing for a numpy integer, so a pandas index
        # value comes back as 'no such spot'. Plain int first.
        try:
            spot_id = int(spot_id)
        except (TypeError, ValueError):
            raise ValueError('%r is not a spot id' % (spot_id,))
        got = self._spot_row(spot_id)
        if got is None:
            raise KeyError('no spot %s in this run' % spot_id)

        gene = str(self.gene_panel[got[0]])
        neighbours, terms, misread, score, prob = self._spot_scores(got)

        rows = []
        for i, internal in enumerate(neighbours[:-1]):
            cls_idx = self._query("SELECT assigned_class_idx FROM cells WHERE %s = ?" % self._cell_key,
                                  (internal,))
            rows.append({
                'cell': self.to_external(internal),
                'class': str(self.class_names[cls_idx[0]]) if cls_idx else None,
                'Gaussian fit': float(terms['mvn_loglik'][i]),
                'class expression': float(terms['attention'][i]),
                'cell scale': float(terms['cell_inefficiency'][i]),
                'cell-gene scale': float(terms['expr_fluct'][i]),
                'gene efficiency': float(terms['gene_inefficiency'][i]),
                'bonus': float(terms['bonus'][i]),
                'sum': float(score[i]),
                'prob': float(prob[i]),
            })
        rows.append({'cell': 'background', 'misread': misread,
                     'sum': misread, 'prob': float(prob[-1])})

        out = {
            'spot': int(spot_id),
            'gene': gene,
            'position': self._spot_position(got),
            'candidates': rows,
            'assigned_to': rows[int(np.argmax(prob))]['cell'],
        }
        out['narrative'] = narrate_spot(out)
        return out

    def cell_counts(self, label, gene=None):
        """How many reads a cell holds. Soft, and the answer says so."""
        c = self._cell(label)
        counts = c['gene_count']
        note = ('these are soft counts: each spot contributes its probability of '
                'belonging to this cell, so they are estimates and not whole numbers')
        if gene is not None:
            hit = np.where(self.gene_panel == gene)[0]
            if not len(hit):
                raise ValueError('no gene %r in the panel' % gene)
            return {'cell': int(label), 'gene': gene,
                    'counts': float(counts[hit[0]]), 'counts_are': note}
        order = np.argsort(-counts, kind='stable')
        return {
            'cell': int(label),
            'total_counts': float(counts.sum()),
            'per_gene': [{'gene': str(self.gene_panel[g]), 'counts': float(counts[g])}
                         for g in order if counts[g] > _COUNT_TOL],
            'counts_are': note,
        }

    def spots_in_cell(self, label, gene=None):
        """How many spots physically sit inside a cell's segmentation mask.

        Hard containment, no probabilities. A spot can sit outside every cell and
        still be assigned to one, so this is a different number from `cell_counts`.
        """
        if not self._has_containment():
            raise NotImplementedError(
                'this run has no containment data. geneData gained the inside_cell '
                'column after it was produced, so the run has to be repeated to '
                'answer this. cell_counts works on any run.')
        label = int(label)
        self.to_internal(label)                  # just to raise if it is not a cell

        import pyarrow.compute as pc

        per_gene = {}
        for tb, gene_of in self._arrow_spots(['inside_cell', 'gene_id']):
            hits = tb.filter(pc.equal(tb['inside_cell'], label))
            for gid, n in zip(*np.unique(hits['gene_id'].to_numpy(), return_counts=True)):
                per_gene[gene_of[int(gid)]] = per_gene.get(gene_of[int(gid)], 0) + int(n)

        note = ('hard count: spots whose pixel falls inside this cell\'s segmentation '
                'mask, no probabilities involved. The model may have assigned some of '
                'them elsewhere, and cell_counts may include spots from outside the mask')
        if gene is not None:
            if gene not in set(gene_of.values()):
                raise ValueError('no gene %r in the panel' % gene)
            return {'cell': label, 'gene': gene, 'spots': per_gene.get(gene, 0),
                    'spots_are': note}
        return {
            'cell': label,
            'total_spots': sum(per_gene.values()),
            # ties by name, so two readers of the same run list them the same way
            'per_gene': [{'gene': g, 'spots': n}
                         for g, n in sorted(per_gene.items(), key=lambda kv: (-kv[1], kv[0]))],
            'spots_are': note,
        }

    def spots_of_cell(self, label, min_prob=None, gene=None):
        """Which spots belong to a cell, under one of two definitions.

        With `min_prob` unset: the spots whose most likely parent is this cell, the
        argmax. With `min_prob` set: every spot with a probability on this cell above
        it, which is what cellData.spot_id holds at 0.0001. The two are different
        lists, and the answer says which one it is. Either way every spot comes with
        its probability, sorted highest first. `gene` keeps only that gene's spots,
        so the number of spots is computed by the tool, not by the agent.
        """
        label = int(label)
        self.to_internal(label)                  # raises if it is not a cell
        if gene is not None:
            self._gene_index(gene)               # raises if it is not a gene

        # the entries that name this cell as a candidate, from the index, in the same
        # order the spot files list them, so the output is what a full scan gave
        idx = self._cell_index()
        spot_row, prob, first = idx.entries_of(label)
        gene_of = idx.gene_of
        rows = []
        for r, p, f in zip(spot_row.tolist(), prob.tolist(), first.tolist()):
            g = gene_of[int(idx.gene_id[r])]
            if gene is not None and g != gene:
                continue
            keep = f if min_prob is None else p > min_prob
            if keep:
                rows.append({'spot': int(idx.spot_id[r]), 'gene': g, 'prob': float(p)})
        rows.sort(key=lambda r: -r['prob'])

        if min_prob is None:
            note = ('spots whose most likely parent is this cell. That is the argmax, '
                    'not a hard assignment: the lowest probability here can be well '
                    'under 0.5, and the cell also draws counts from spots whose most '
                    'likely parent is another cell')
        else:
            note = ('every spot with probability above %g on this cell. Their '
                    'probabilities add up to the cell\'s soft counts, so this is the '
                    'decomposition of cell_counts. Probabilities are rounded to 3 '
                    'decimals in the viewer files, so anything under 0.0005 is absent'
                    % min_prob)
        return {
            'cell': label,
            'gene': gene,
            'definition': 'most likely parent' if min_prob is None else 'prob > %g' % min_prob,
            'n_spots': len(rows),
            'sum_of_probs': float(sum(r['prob'] for r in rows)),
            'spots': rows,
            'spots_are': note,
        }

    # ----------------------------------------------------- the whole run

    def _all_cells(self):
        """Every cell's class probabilities, gene counts and assigned class, read
        once and kept. Rows are in internal order, the background row left out.
        A few tens of MB on a big run, which is fine for the questions these
        answer (cells per class, where is a gene, which cells match)."""
        if not hasattr(self, '_cells_scan'):
            with self._lock:
                rows = self._con.execute(
                    'SELECT %s, class_prob, gene_count, assigned_class_idx FROM cells '
                    'WHERE %s > 0 ORDER BY %s' % ((self._cell_key,) * 3)).fetchall()
            self._cells_scan = {
                'internal': np.array([r[0] for r in rows], dtype=np.int64),
                'class_prob': np.stack([np.frombuffer(r[1], dtype=np.float32) for r in rows]),
                'gene_count': np.stack([np.frombuffer(r[2], dtype=np.float32) for r in rows]),
                'assigned': np.array([r[3] for r in rows], dtype=np.int64),
            }
        return self._cells_scan

    def _positions(self):
        """Fitted cell centroid (X, Y, Z) per segmentation label from the viewer's cells table,
        or None when the viewer files are not in this run."""
        if not hasattr(self, '_pos'):
            try:
                tb = self._arrow_cells()
                self._pos = {int(c): (float(x), float(y), float(z)) for c, x, y, z in zip(
                    tb['cell_id'].to_pylist(), tb['X'].to_pylist(),
                    tb['Y'].to_pylist(), tb['Z'].to_pylist())}
            except Exception:
                self._pos = None
        return self._pos

    def _gene_index(self, name):
        hit = np.where(self.gene_panel == name)[0]
        if not len(hit):
            raise ValueError('no gene %r in this run' % name)
        return int(hit[0])

    def _class_index(self, name):
        hit = np.where(self.class_names == name)[0]
        if not len(hit):
            raise ValueError('no class %r in this run' % name)
        return int(hit[0])

    def gene(self, name):
        """One gene across the run: its efficiency, misread density, how many spots
        it has and where they went, class by class and cell by cell."""
        g = self._gene_index(name)
        cfg = self.config or {}
        eta = float(np.asarray(self._meta('eta_bar', parse=True), dtype=np.float32)[g])
        scan = self._all_cells()
        col = scan['gene_count'][:, g]
        soft_per_class = scan['class_prob'].T @ col
        hard_per_class = np.bincount(scan['assigned'], weights=col, minlength=self.nK)
        order = np.argsort(-soft_per_class, kind='stable')
        top = np.argsort(-col, kind='stable')[:10]
        totals = self._meta('gene_total_spots', parse=True, default=None)
        misread = self._meta('hard_misread_counts', parse=True, default=None)
        return {
            'gene': name,
            'eta': eta,
            'eta_is': 'eta_bar, the posterior mean of eta, the gene\'s detection efficiency '
                      'relative to the Inefficiency setting; the reference '
                      'expression of every class is multiplied by it',
            'eta_times_Inefficiency': eta * cfg['Inefficiency'] if cfg.get('Inefficiency') is not None else None,
            'eta_times_Inefficiency_is': 'eta_bar times the Inefficiency setting, as Genes.inefficiency '
                               'defines it; None on runs without the settings record',
            'misread_density': self._meta('rho_bar', parse=True, default={}).get(name),
            'misread_density_prior': self._meta('misread_density', parse=True, default={}).get(name),
            'misread_density_is': 'rho_bar, the rate of misread spots per unit volume the '
                                  'model learned for this gene; the prior is where it started',
            'total_spots': int(totals[g]) if totals else None,
            'spots_called_misread': int(misread[g]) if misread else None,
            'counts_in_cells': float(col.sum()),
            'counts_are': 'soft, weighted by the spot assignment probabilities',
            'counts_per_class': [{'class': str(self.class_names[k]),
                                  'soft': float(soft_per_class[k]),
                                  'in_cells_called': float(hard_per_class[k])}
                                 for k in order if soft_per_class[k] > _COUNT_TOL],
            'counts_per_class_are': 'soft: every cell weighted by its probability of the '
                                    'class; in_cells_called: summed over the cells whose '
                                    'most probable class it is',
            'top_cells': [{'cell': self.to_external(int(scan['internal'][i])),
                           'counts': float(col[i]),
                           'class': str(self.class_names[scan['assigned'][i]])}
                          for i in top if col[i] > _COUNT_TOL],
        }

    def theta(self, label):
        """The cell scale factor of one cell, overall and under each class."""
        c = self._cell(label)
        order = np.argsort(-c['class_prob'], kind='stable')
        return {
            'cell': int(label),
            'theta': c['theta'],
            'theta_is': 'the scale factor averaged over the class probabilities, the '
                        'value the outputs report',
            'theta_bar_per_class': [{'class': str(self.class_names[k]),
                                     'theta_bar': float(c['theta_bar'][k]),
                                     'prob': float(c['class_prob'][k])}
                                    for k in order[:5]],
            'theta_bar_is': 'the posterior mean of theta under each class: the factor '
                            'scaling that class\'s expected counts to the cell\'s total. '
                            'Gamma(rTheta, rTheta) prior, mean 1',
            'rTheta': (self.config or {}).get('rTheta'),
        }

    def gamma(self, label, gene=None):
        """The cell-gene scale factors of one cell under its assigned class."""
        c = self._cell(label)
        gam = c['gamma_assigned']
        cls = str(self.class_names[c['assigned_class_idx']])
        out = {
            'cell': int(label),
            'class': cls,
            'gamma_is': 'gamma_bar, the posterior mean of the per cell, per gene scale that '
                        'absorbs overdispersion, under the assigned class only. '
                        'diagnostics.db does not keep the (cell, gene, class) array, so '
                        'gamma under another class is not available',
            'expected_is': 'the count the assigned class predicts for the gene in this cell, '
                           'after Inefficiency, eta and theta. gamma is not counts over '
                           'expected: it is (rSpot + counts) over (rSpot + expected), which '
                           'pulls it toward 1, so say how far off a gene is from counts '
                           'against expected, not from gamma',
        }
        # gamma_bar = (rSpot + counts) / (rSpot + expected), see VarBayes.gamma_upd.
        # Turned round it gives back the expected count the fit used
        r_spot = float(self._meta('rSpot'))
        expected = np.maximum((r_spot + c['gene_count']) / gam - r_spot, 0)
        if gene is not None:
            g = self._gene_index(gene)
            out['gene'] = gene
            out['gamma'] = float(gam[g])
            out['counts'] = float(c['gene_count'][g])
            out['expected'] = float(expected[g])
        else:
            out['gamma'] = [{'gene': str(self.gene_panel[g]), 'gamma': float(gam[g]),
                             'counts': float(c['gene_count'][g]),
                             'expected': float(expected[g])}
                            for g in range(self.nG)]
        return out

    def spot(self, spot_id):
        """One spot: its gene, position, plane, and the cells it may belong to
        with their probabilities. explain_spot gives the terms behind them."""
        try:
            spot_id = int(spot_id)
        except (TypeError, ValueError):
            raise ValueError('%r is not a spot id' % (spot_id,))
        got = self._spot_row(spot_id)
        if got is None:
            raise KeyError('no spot %s in this run' % spot_id)
        neighbours, _, _, _, prob = self._spot_scores(got)
        cands = []
        for i, internal in enumerate(neighbours[:-1]):
            k = self._query("SELECT assigned_class_idx FROM cells WHERE %s = ?" % self._cell_key,
                            (internal,))
            cands.append({'cell': self.to_external(internal),
                          'class': str(self.class_names[k[0]]) if k else None,
                          'prob': float(prob[i])})
        cands.append({'cell': 'background', 'prob': float(prob[-1])})
        best = int(np.argmax(prob))
        return {
            'spot': spot_id,
            'gene': str(self.gene_panel[got[0]]),
            'position': self._spot_position(got),
            'assigned_to': cands[best]['cell'],
            'prob': float(prob[best]),
            'candidates': cands,
            'candidates_are': 'the nearest cells the spot was scored against, plus the '
                              'background, with the probability of each',
        }

    def neighbours(self, label):
        """The cells whose classes enter this cell's spatial term, with their
        classes and, when the viewer files are there, their distances."""
        label = int(label)
        nbrs = self._neighbours(label)
        if nbrs is None:
            raise NotImplementedError('this run was written before diagnostics.db kept the '
                                      'spatial neighbours; rerun with the current pciSeq')
        scan = self._all_cells()
        pos = self._positions()
        me = pos.get(label) if pos else None
        rows = []
        for n in nbrs:
            i = int(np.searchsorted(scan['internal'], self.to_internal(n)))
            k = int(scan['assigned'][i])
            r = {'cell': n, 'class': str(self.class_names[k]),
                 'prob': float(scan['class_prob'][i, k])}
            if me and n in pos:
                x, y, z = pos[n]
                r['distance_xy'] = float(np.hypot(x - me[0], y - me[1]))
                r['plane_offset'] = float(z - me[2])
            rows.append(r)
        c = self._cell(label)
        return {
            'cell': label,
            'class': str(self.class_names[c['assigned_class_idx']]),
            'n': len(rows),
            'neighbours': rows,
            'neighbours_are': 'the cells the spatial (mrf) term listens to, nearest first. '
                              'distance_xy is between fitted cell centroids in pixels of the xy plane, '
                              'plane_offset in planes',
            'mrf_beta': (self.config or {}).get('mrf_beta'),
        }

    def class_counts(self, min_counts=None):
        """How many cells each class has, hard and soft."""
        scan = self._all_cells()
        totals = scan['gene_count'].sum(axis=1)
        keep = totals >= min_counts if min_counts is not None else np.ones(len(totals), bool)
        hard = np.bincount(scan['assigned'][keep], minlength=self.nK)
        soft = scan['class_prob'][keep].sum(axis=0)
        rows = [{'class': str(name), 'cells': int(hard[k]), 'soft': float(soft[k])}
                for k, name in enumerate(self.class_names)]
        zero = [r for r in rows if r['class'] == 'Zero']
        rest = sorted([r for r in rows if r['class'] != 'Zero'], key=lambda r: -r['cells'])
        return {
            'n_cells': int(keep.sum()),
            'min_counts': min_counts,
            'classes': zero + rest,
            'cells_is': 'hard: the number of cells whose most probable class it is',
            'soft_is': 'the class probability summed over the cells, the expected '
                       'number of cells of the class',
            'min_counts_is': 'cells with fewer soft counts in total are left out',
        }

    def class_prior(self):
        """The weight of each class at the end of the run.

        How likely a cell is to be of the class before any of its genes is looked
        at. Read from the log_prior the run saved, which is the prior term the class
        score adds, and not rebuilt from the cell counts: the counts give the mean
        of the weight, while the score uses the mean of its log, which is lower for
        a class with only a few cells.
        """
        log_prior = self._meta('log_prior', parse=True, default=None)
        if log_prior is None:
            raise KeyError('this run did not save log_prior, so its class weights are not known')
        rows = [{'class': str(name), 'weight': float(np.exp(log_prior[k])),
                 'log_weight': float(log_prior[k])}
                for k, name in enumerate(self.class_names)]
        zero = [r for r in rows if r['class'] == 'Zero']
        rest = sorted([r for r in rows if r['class'] != 'Zero'], key=lambda r: -r['weight'])
        return {
            'cell_type_prior': (self.config or {}).get('cell_type_prior'),
            'classes': zero + rest,
            'weight_is': 'the prior probability of the class at the end of the run, the one '
                         'the class score of a cell uses: exp of log_weight',
            'cell_type_prior_is': 'weighted: the weights were learned from the cells and '
                                  'refitted at every iteration, Zero apart, whose weight is '
                                  'fixed. uniform: all of them were fixed by cell_type_weights',
            'note': 'a cell that was pinned partway through the run kept the weights of that '
                    'iteration, so its own prior term can differ from these',
        }

    def find_cells(self, class_name=None, plane=None, min_counts=None,
                   top_two_within=None, n=50):
        """The cells matching a few filters: assigned class, plane of the fitted cell centroid,
        total counts, and how close the top two classes are."""
        scan = self._all_cells()
        mask = np.ones(len(scan['internal']), bool)
        if class_name is not None:
            mask &= scan['assigned'] == self._class_index(class_name)
        totals = scan['gene_count'].sum(axis=1)
        if min_counts is not None:
            mask &= totals >= min_counts
        top2 = -np.sort(-scan['class_prob'], axis=1)[:, :2]
        margin = top2[:, 0] - top2[:, 1]
        if top_two_within is not None:
            mask &= margin <= top_two_within
        labels = [self.to_external(int(i)) for i in scan['internal']]
        if plane is not None:
            pos = self._positions()
            if pos is None:
                raise NotImplementedError('the viewer files are not in this run, so the '
                                          'cells have no positions to filter by plane')
            if self._plane_of(0.0) is None:
                raise NotImplementedError('this run does not record its voxel size, so a '
                                          'fitted cell centroid cannot be put on a plane; rerun with the '
                                          'current pciSeq to record it')
            # Z in the viewer files is the anisotropy scaled z, not a plane index
            planes = np.array([self._plane_of(pos[lab][2]) if lab in pos else -1
                               for lab in labels])
            mask &= planes == int(plane)
        idx = np.where(mask)[0]
        idx = idx[np.argsort(-top2[idx, 0], kind='stable')]
        runner = np.argsort(-scan['class_prob'], axis=1)[:, 1]
        rows = [{'cell': labels[i],
                 'class': str(self.class_names[scan['assigned'][i]]),
                 'prob': float(top2[i, 0]),
                 'runner_up': str(self.class_names[runner[i]]),
                 'margin': float(margin[i]),
                 'total_counts': float(totals[i])}
                for i in idx[:n]]
        return {
            'n_matching': int(len(idx)),
            'shown': len(rows),
            'filters': {'class_name': class_name, 'plane': plane, 'min_counts': min_counts,
                        'top_two_within': top_two_within},
            'cells': rows,
            'cells_are': 'matched on the most probable class; margin is the probability of '
                         'the assigned class minus the runner up; plane is the plane the '
                         'fitted cell centroid falls on, rounded down as for the spots; sorted by prob, '
                         'the first n shown',
        }

    def metadata(self, key=None):
        """The metadata table of diagnostics.db: the keys, or one value."""
        with self._lock:
            rows = self._con.execute('SELECT key, value FROM metadata').fetchall()

        def parsed(v):
            try:
                return json.loads(v)
            except (TypeError, ValueError):
                return v

        def kind(v):
            v = parsed(v)
            if isinstance(v, list):
                inner = ' of lists of %d' % len(v[0]) if v and isinstance(v[0], list) else ''
                return 'list of %d%s' % (len(v), inner)
            if isinstance(v, dict):
                return 'dict of %d keys' % len(v)
            return 'number' if isinstance(v, (int, float)) else 'text'

        if key is None:
            return {'keys': [{'key': k, 'kind': kind(v)} for k, v in rows],
                    'note': 'call again with a key for its value. eta_bar, gene_total_spots '
                            'and hard_misread_counts follow gene_panel; log_prior and '
                            'class_names go together; mean_gene_reads_per_class and '
                            'sc_mean_expression are gene by class'}
        for k, v in rows:
            if k == key:
                return {'key': key, 'value': parsed(v)}
        raise KeyError('no metadata key %r, call metadata() for the list' % key)

    def cell_row(self, label):
        """The cellData row of one cell, value for value.

        Everything the tsv has except the two drawing columns. Read straight out of
        cellData.tsv when it is there, which is both exact and fast, a targeted grep
        on a 250 MB file takes about 30 ms. Without the tsv the row is rebuilt from
        the viewer files, which matches on every column except the boundary cases of
        spot_id; see `_exact_spot_ids`.
        """
        label = int(label)
        self.to_internal(label)

        row = self._tsv_row('cellData.tsv', r'^%d\t' % label)
        if row is not None:
            # X0, Y0, Z0 is the segmentation centroid, runs before October 2026
            # do not have it, hence the `if k in row`
            out = {k: row[k] for k in ('Cell_Num', 'X', 'Y', 'Z', 'X0', 'Y0', 'Z0',
                                       'Genenames', 'CellGeneCount', 'spot_id',
                                       'ClassName', 'Prob') if k in row}
            out['source'] = 'cellData.tsv'
        else:
            import pyarrow.compute as pc
            tb = self._arrow_cells()
            hit = tb.filter(pc.equal(tb['cell_id'], label))
            if hit.num_rows == 0:
                raise KeyError('no cell %s in this run' % label)
            r = hit.to_pylist()[0]
            genes = r['gene_names']
            per_gene = self._exact_spot_ids(label)
            out = {
                'Cell_Num': label,
                'X': _tsv(r['X']), 'Y': _tsv(r['Y']), 'Z': _tsv(r['Z']),
                'Genenames': genes,
                'CellGeneCount': [_tsv(v) for v in r['gene_counts']],
                'spot_id': [per_gene.get(g, []) for g in genes],
                'ClassName': r['class_name'],
                'Prob': [_tsv(v) for v in r['prob']],
                'source': 'rebuilt from the viewer files, cellData.tsv is not in this '
                          'run. spot_id can miss a few spots sitting within 0.0005 of '
                          'the 0.0001 cut-off',
            }
        out['counts_are'] = ('soft, weighted by the spot assignment probabilities. '
                             'spot_id lists every spot with probability above 0.0001 on '
                             'this cell, lined up with Genenames, and their probabilities '
                             'add up to CellGeneCount')
        # the plane the centroid is on, given so the agent never works it out: Z is
        # the anisotropy scaled z, and a model read Z = 9.256 of espio cell 8 as
        # plane 9 when it is plane 3. From the full z in diagnostics.db, Z above is
        # cut to 3 decimals and could tip a centroid on a plane boundary
        got = self._query('SELECT z FROM cells WHERE %s = ?' % self._cell_key,
                          (self.to_internal(label),))
        plane = self._plane_of(got[0]) if got and got[0] is not None else None
        if plane is not None:
            out['plane'] = plane
            out['Z_is'] = ('the anisotropy scaled z the model works in, not a plane; plane '
                           'is the plane the fitted cell centroid sits on')
        else:
            out['Z_is'] = ('the anisotropy scaled z the model works in, not a plane. This '
                           'run carries no voxel_size, so the plane cannot be given')
        return out

    def spot_row(self, spot_id):
        """The geneData row of one spot, value for value.

        Read out of geneData.tsv when it is there, otherwise rebuilt from the viewer
        files; `source` says which. Both give the same numbers, the viewer files are
        written from the same frame.

        neighbour_prob is rounded to 3 decimals here because that is how geneData
        stores it (summary.py:111), so a candidate holding 0.0004 shows up as 0.0 and
        you cannot tell it apart from one holding nothing at all. That is fine for
        reporting what the file says, and it is why cell_row's spot_id list can
        contain spots whose probability reads 0.0: the cut-off there is 0.0001,
        applied before the rounding. When the small probabilities matter use
        explain_spot, which recomputes them from diagnostics.db at full precision.
        """
        try:
            spot_id = int(spot_id)
        except (TypeError, ValueError):
            raise ValueError('%r is not a spot id' % (spot_id,))
        # spot_id is the third column of geneData
        row = self._tsv_row('geneData.tsv', r'^[^\t]*\t[^\t]*\t%d\t' % spot_id)
        if row is not None:
            row['source'] = 'geneData.tsv'
            return row

        import pyarrow.compute as pc
        for tb, gene_of in self._arrow_spots(None):
            hit = tb.filter(pc.equal(tb['spot_id'], spot_id))
            if hit.num_rows:
                r = hit.to_pylist()[0]
                out = {
                    'gene_name': gene_of[int(r['gene_id'])],
                    'gene_id': int(r['gene_id']),
                    'spot_id': spot_id,
                    'x': _tsv(r['x']), 'y': _tsv(r['y']), 'z': _tsv(r['z']),
                    'plane_id': int(r['plane_id']),
                    'neighbour': int(r['neighbour_array'][0]),
                    'neighbour_array': [int(c) for c in r['neighbour_array']],
                    'neighbour_prob': [_tsv(v) for v in r['neighbour_prob']],
                    'omp_score': _tsv(r['omp_score']),
                    'omp_intensity': _tsv(r['omp_intensity']),
                    'is_hard_misread': int(r['is_hard_misread']),
                    'source': 'the viewer files, geneData.tsv is not in this run',
                }
                if 'inside_cell' in r:
                    out['inside_cell'] = int(r['inside_cell'])
                return out
        raise KeyError('no spot %s in this run' % spot_id)

    # ------------------------------------------------------------ internals

    def _mean_expression_adj(self):
        """The (nG, nK) matrix the likelihood scores against, from the metadata of a
        format 2 file."""
        if self._mea is None:
            self._mea = np.asarray(self._meta('mean_expression_adj', parse=True), dtype=np.float32)
        return self._mea

    def _spot_row(self, spot_id):
        """One row of the spots table as format 1 laid it out, whatever the file's
        version: (gene_idx, x, y, z, the nN candidate ids with the background 0 last,
        and the six terms as float32 arrays of nN with the background slot zero).
        None when the spot is not in the run."""
        has = self._spot_cols
        cols = ['gene_idx', 'x', 'y', 'z', 'neighbor_cell_ids', 'mvn_loglik', 'attention',
                'expr_fluct', 'cell_inefficiency']
        cols += ['gene_inefficiency'] if 'gene_inefficiency' in has else []
        cols += ['bonus'] if 'bonus' in has else []
        got = self._query("SELECT %s FROM spots WHERE spot_id = ?" % ', '.join(cols), (spot_id,))
        if got is None:
            return None
        got = dict(zip(cols, got))
        ids = got['neighbor_cell_ids']
        if isinstance(ids, (bytes, memoryview)):
            # format 2: the cells only, int32
            ids = np.frombuffer(ids, dtype=np.int32).tolist() + [0]
        else:
            ids = json.loads(ids)
        n = len(ids)
        f32 = lambda b: np.frombuffer(b, dtype=np.float32)
        pad = lambda a: a if len(a) == n else np.append(a, np.zeros(n - len(a), dtype=np.float32))
        terms = [pad(f32(got[k])) for k in ('mvn_loglik', 'attention', 'expr_fluct', 'cell_inefficiency')]
        if 'gene_inefficiency' in got:
            terms.append(pad(f32(got['gene_inefficiency'])))
        else:
            le = np.float32(self._meta('logeta_bar', parse=True)[got['gene_idx']])
            terms.append(np.append(np.full(n - 1, le, dtype=np.float32), np.float32(0)))
        terms.append(pad(f32(got['bonus'])) if 'bonus' in got else np.zeros(n, dtype=np.float32))
        return (got['gene_idx'], got['x'], got['y'], got['z'], ids, *terms)

    def _spot_scores(self, got):
        """The candidate ids, the six terms, the misread density, the scores and the
        softmax for one spots-table row as _spot_row lays it out."""
        gene = str(self.gene_panel[got[0]])
        neighbours = got[4]
        terms = {name: got[i]
                 for i, name in enumerate(
                     ['mvn_loglik', 'attention', 'expr_fluct',
                      'cell_inefficiency', 'gene_inefficiency', 'bonus'], start=5)}
        misread = float(self._meta('log_rho_bar', parse=True)[gene])
        # the last slot of every array is the background, which scores the misread
        # density instead of the six terms
        score = sum(v[:-1] for v in terms.values())
        score = np.append(score, misread)
        prob = np.exp(score - score.max())
        prob = prob / prob.sum()
        return neighbours, terms, misread, score, prob

    def _exact_spot_ids(self, label, tol=0.0001):
        """cellData.spot_id for one cell: per gene, the spots with an unrounded
        probability above tol on the cell.

        The arrow files round probabilities to 3 decimals, so they cannot answer this
        on their own, but they do list every candidate of every spot whatever its
        probability. So: candidates from the arrow files, then the exact probability
        of each from the db, which is indexed by spot id.

        This is the fallback, only used when cellData.tsv is missing, and it is not
        quite exact. The recomputation sits within about 5e-4 of what the model held,
        which is nothing anywhere except at the 0.0001 cut-off, where it flips a spot
        in or out. Measured on espio cell 18223: 506 spots against the 542 cellData
        lists, the 36 missing ones all recomputing to just under 1e-4. Read the tsv
        when it is there.
        """
        internal = self.to_internal(label)
        idx = self._cell_index()
        spot_row, _, _ = idx.entries_of(int(label))
        candidates = [int(idx.spot_id[r]) for r in spot_row.tolist()]

        per_gene = {}
        for sid in candidates:
            got = self._spot_row(sid)
            neighbours, _, _, _, prob = self._spot_scores(got)
            p = float(prob[neighbours.index(internal)])
            if p > tol:
                per_gene.setdefault(str(self.gene_panel[got[0]]), []).append(sid)
        return per_gene

    def _tsv_row(self, name, pattern):
        """One row of a tsv, found by grep, parsed with the header.

        None when the file is not there. grep rather than pandas on purpose: these
        files run to hundreds of MB and reading one row out of cellData.tsv this way
        takes about 30 ms, against seconds to load the lot.
        """
        import ast
        import subprocess
        path = self.path.parent.parent.parent / 'tsv' / name
        if not path.is_file():
            return None
        hit = subprocess.run(['grep', '-m1', '-P', pattern, str(path)],
                             capture_output=True, text=True)
        if hit.returncode or not hit.stdout.strip():
            return None
        with open(path) as fid:
            header = fid.readline().rstrip('\n').split('\t')

        out = {}
        for col, raw in zip(header, hit.stdout.rstrip('\n').split('\t')):
            if raw.startswith('[') or raw.startswith('{'):
                out[col] = ast.literal_eval(raw)
            else:
                try:
                    out[col] = int(raw)
                except ValueError:
                    try:
                        out[col] = float(raw)
                    except ValueError:
                        out[col] = raw
        return out

    def _arrow_cells(self):
        """The arrow cells table, read once. It is cellData minus spot_id and the
        drawing columns, at the tsv's precision."""
        if not hasattr(self, '_cells_tb'):
            import pyarrow as pa
            import pyarrow.feather as feather
            d = self.path.parent.parent / 'arrow_cells'
            self._cells_tb = pa.concat_tables(
                [feather.read_table(f) for f in sorted(d.glob('*.feather'))])
        return self._cells_tb

    def _cell_index(self):
        """Every spot's candidate cells, grouped by cell, built once per run.

        spots_of_cell and the cellData fallback ask "which spots name this cell as a
        candidate". Scanning the 3 million spots of a big run for that took close to
        a minute per question; this reads the spot files once, a few seconds, and
        then each question touches only that cell's own entries.

        Held: per spot, its id and gene id; per candidate entry, grouped by cell in
        the order the files list them, the spot it belongs to, its probability and
        whether it is the spot's first, most likely, candidate. About 9 bytes per
        entry, around 300 MB on espio, only for a run that is asked these questions.
        """
        if getattr(self, '_cidx', None) is None:
            with self._lock:
                if getattr(self, '_cidx', None) is None:
                    self._cidx = _CellIndex(self._arrow_spots(
                        ['spot_id', 'gene_id', 'neighbour_array', 'neighbour_prob']))
        return self._cidx

    def _arrow_spots(self, columns):
        """Yield each arrow spots shard as a table, plus the gene id to name map.

        The shards carry segmentation labels, they went through the same translation
        as geneData, so a label compares directly against neighbour_array.
        """
        import pyarrow.feather as feather
        shard_dir = self.path.parent.parent / 'arrow_spots'
        gene_of = {int(k): v for k, v in
                   json.loads((shard_dir / 'gene_dict.json').read_text()).items()}
        for f in sorted(shard_dir.glob('*.feather')):
            yield feather.read_table(f, columns=columns), gene_of

    def _background(self, channel=None, mbtiles=None):
        """Which background image to draw on, as (path, its name, every name there is,
        a note for the answer or None).

        Looked up the way the viewer does it: every .mbtiles in viewer_data is one
        background image, named by the 'name' in its metadata, or by its file name if
        it has none. `channel` picks one by that name or by the file name, ignoring
        case; a part of the name is enough if it fits only one image ('dapi' finds
        'dapi_Espio'). When there are several images and no channel we do not guess,
        we refuse and list them, so the agent asks the user.

        With only one image the name is not checked at all. Name and description are
        free text, often empty, so they cannot tell us whether it is DAPI or GCaMP;
        refusing 'show me it on the DAPI' because the file is called 'Espio' would
        be wrong. We draw on it and say so in the note.
        """
        if mbtiles:
            path = Path(mbtiles).expanduser()
            return path, _channel_name(path), [_channel_name(path)], None
        viewer_data = self.path.parent.parent
        found = sorted(viewer_data.glob('*.mbtiles'))
        if not found:
            raise FileNotFoundError(
                'no .mbtiles in %s, so there is no image to draw on. The viewer tiles are '
                'made by pciSeq.stage_image; copy the file into viewer_data or pass its '
                'path as mbtiles' % viewer_data)
        names = [_channel_name(f) for f in found]
        listing = ', '.join('%s (%s)' % (n, f.name) for n, f in zip(names, found))

        if len(found) == 1:
            note = None
            if channel is not None:
                note = ('asked for %r, but this run has only one background image, '
                        'called %r, so that is the one drawn. Which stain it is cannot be '
                        'told from the file, only whoever made it knows'
                        % (channel, names[0]))
            return found[0], names[0], names, note
        if channel is None:
            raise ValueError('this run has %d background images: %s. Ask the user which '
                             'one they want and pass it as channel' % (len(found), listing))

        want = str(channel).strip().lower()
        keys = [(n.lower(), f.stem.lower()) for n, f in zip(names, found)]
        hits = [i for i, k in enumerate(keys) if want in k]
        if not hits:
            hits = [i for i, k in enumerate(keys) if any(want in s for s in k)]
        if len(hits) != 1:
            raise ValueError('%s background image matching %r. This run has: %s'
                             % ('no' if not hits else 'more than one', channel, listing))
        return found[hits[0]], names[hits[0]], names, None

    def _outlines(self, plane=None):
        """The cell outlines, from the viewer's per plane boundary files.

        With a plane, every outline on that plane as a table. Without one, the whole
        folder, which is what you need to find which planes a cell is on.
        """
        import pyarrow as pa
        import pyarrow.feather as feather
        d = self.path.parent.parent / 'arrow_boundaries'
        if plane is not None:
            f = d / ('boundaries_plane_%02d.feather' % plane)
            if not f.exists():
                raise ValueError('no outlines for plane %s in this run' % plane)
            return feather.read_table(f)
        files = sorted(d.glob('boundaries_plane_*.feather'))
        if not files:
            raise FileNotFoundError('no arrow_boundaries in %s' % d.parent)
        return pa.concat_tables([feather.read_table(f) for f in files])

    def _neighbours(self, label):
        """The cells the mrf term of this cell listens to, as segmentation labels,
        nearest first. None on runs written before diagnostics.db kept them."""
        row = self.to_internal(label)
        if not self._has_nbrs:
            return None
        got = self._query('SELECT neighbours FROM cells WHERE %s = ?' % self._cell_key,
                          (row,))
        return [self.to_external(r) for r in np.frombuffer(got[0], dtype=np.int32)]

    def cell_image(self, label, context=False, plane=None, width=1200, channel=None,
                   neighbours=False, mbtiles=None):
        """A picture of one cell on the background image of the run.

        The background is stitched from the viewer's tile pyramid with read_tiles, so
        it is a close copy of the image, not the raw pixels (the tiles are jpeg). Two
        kinds of picture:

        * close-up (context=False): the cell outlined in red and every other cell on
          that plane in blue, in a 3:2 window about four and a half times the size of
          the cell. With neighbours=True only the cells the mrf term listens to are
          outlined, the picture for 'why did its neighbours make it this class'. Those
          on another plane have no outline here and the answer lists them.
        * context (context=True): the whole plane trimmed to 3:2 with a ring round the
          cell, to show where in the tissue it sits.

        The plane is the one given, else the plane of the fitted cell centroid when the run
        carries voxel_size, else the plane where the cell outline is biggest. Returns
        the image (PIL) and a dict of what was drawn. The MCP tool also takes
        save_as; here just call .save on the image.

        channel picks the background image when the run has more than one (DAPI,
        GCaMP, ...); with several and no channel it refuses and lists them. With one
        image channel is not checked, the one image is used.
        """
        import pyarrow.compute as pc
        from PIL import ImageDraw, ImageOps
        from ..tiling.read_tiles import read_tiles

        label = int(label)
        if self.to_internal(label) == 0:
            raise ValueError('cell 0 is the background pseudocell, not a cell')
        mbtiles, channel, channels, note = self._background(channel, mbtiles)

        allb = self._outlines()
        mine = allb.filter(pc.equal(allb['label'], label)).to_pylist()
        if not mine:
            raise KeyError('cell %s has no outline in arrow_boundaries' % label)
        area = {r['plane_id']: _polygon_area(r['x_list'], r['y_list']) for r in mine}

        cells = self._arrow_cells()
        hit = cells.filter(pc.equal(cells['cell_id'], label)).to_pylist()
        cx, cy, cz = (hit[0]['X'], hit[0]['Y'], hit[0]['Z']) if hit else (None, None, None)

        if plane is not None:
            plane, why = int(plane), 'asked for'
        else:
            plane = self._plane_of(cz) if cz is not None else None
            why = 'the plane of the fitted cell centroid'
            if plane not in area:
                # old runs have no voxel_size, so take the biggest slice through the cell
                plane = max(area, key=area.get)
                why = 'the plane where the cell outline is biggest'
        outline = next((r for r in mine if r['plane_id'] == plane), None)

        img_w, img_h, _ = _mbtiles_meta(mbtiles)

        # the centre is the centroid, or the middle of the outline if the cell is not in
        # arrow_cells for some reason
        if cx is None:
            cx, cy = np.mean(mine[0]['x_list']), np.mean(mine[0]['y_list'])
        red, blue = (255, 96, 80), (120, 210, 255)

        if context:
            # whole plane cut down to 3:2, trimming whichever side is too long
            if img_w / img_h > 1.5:
                bw, bh = img_h * 1.5, img_h
            else:
                bw, bh = img_w, img_w / 1.5
            box = ((img_w - bw) / 2, (img_h - bh) / 2, (img_w + bw) / 2, (img_h + bh) / 2)
            im, scale = read_tiles(str(mbtiles), plane=plane, bbox=box, width=width)
            im = ImageOps.autocontrast(im, cutoff=(0.2, 0.05))  # whole plane is dark
            ax, ay = (cx - box[0]) * scale, (cy - box[1]) * scale
            r = 0.022 * max(im.size)
            ImageDraw.Draw(im).ellipse([ax - r, ay - r, ax + r, ay + r], outline=red,
                                       width=max(2, width // 200))
            others = 0
        else:
            # window 4.5 times the cell, 3:2. For cell 18223 on espio that is the
            # 210 x 140 of the docs figure
            ref = outline or max(mine, key=lambda r: area[r['plane_id']])
            ext = max(np.ptp(ref['x_list']), np.ptp(ref['y_list']), 10)
            bh = min(4.5 * ext, img_h)
            bw = min(1.5 * bh, img_w)
            # keep the box inside the image rather than letting read_tiles clip it, so
            # the origin used for the outlines below stays right
            x0 = min(max(cx - bw / 2, 0), img_w - bw)
            y0 = min(max(cy - bh / 2, 0), img_h - bh)
            box = (x0, y0, x0 + bw, y0 + bh)
            im, scale = read_tiles(str(mbtiles), plane=plane, bbox=box, width=width)

            on_plane = self._outlines(plane).to_pylist()
            nbrs = self._neighbours(label) if neighbours else None
            if nbrs is not None:
                near = [r for r in on_plane if r['label'] in nbrs]
            else:
                near = [r for r in on_plane if r['label'] != label
                        and max(r['x_list']) >= box[0] and min(r['x_list']) <= box[2]
                        and max(r['y_list']) >= box[1] and min(r['y_list']) <= box[3]]
            thin = max(2, width // 300)
            _draw_outlines(im, scale, box, near, blue, thin, 0.15)
            if outline:
                _draw_outlines(im, scale, box, [outline], red, thin + 2, 0.3)
            others = len(near)

        if context:
            picture = 'context, the whole plane with the cell ringed'
        elif neighbours and nbrs is not None:
            picture = 'close-up, the cell in red and its mrf neighbours in blue'
        else:
            picture = 'close-up, the cell in red and the other cells on the plane in blue'
        extra = {}
        if neighbours and not context:
            if nbrs is None:
                extra['neighbours_note'] = (
                    'this run was written before diagnostics.db kept the mrf '
                    'neighbours, so every cell in the window is outlined instead. '
                    'Rerunning with the current pciSeq records them')
            else:
                drawn = {r['label'] for r in near}
                extra['neighbours'] = nbrs
                extra['neighbours_not_on_this_plane'] = [n for n in nbrs if n not in drawn]
        info = {
            'cell': label,
            'picture': picture,
            'plane': plane,
            'plane_is': why,
            # in words: as a list [12, 34] a model read it as 'planes 12 and 34'
            'outline_on_planes': '%d to %d' % (min(area), max(area)),
            'cell_has_outline_on_this_plane': outline is not None,
            'centroid_xy': [float(cx), float(cy)],
            'bbox': [round(float(v), 1) for v in box],
            'scale': round(float(scale), 3),
            'other_cells_outlined': others,
            **extra,
            'background': channel,
            'backgrounds_in_this_run': channels,
            **({'background_note': note} if note else {}),
            'mbtiles': str(mbtiles),
            'image_is': 'stitched from the jpeg tile pyramid, a close visual copy of the '
                        'image, not the raw pixels. Fine to look at, not to measure',
        }
        return im, info

    def plane_image(self, plane=None, bbox=None, width=1200, channel=None, mbtiles=None):
        """The background image of one plane, whole or a part of it, with nothing
        drawn on top.

        Stitched from the viewer's tile pyramid with read_tiles, so like cell_image
        it is a close copy of the image, not the raw pixels. The plane defaults to
        the middle one of the stack. bbox is (x0, y0, x1, y1) in image pixels, the
        same coordinates as the cells and spots, and is clamped to the image; leave
        it out for the whole plane, untrimmed. Returns the image (PIL) and a dict of
        what was read. The MCP tool also takes save_as; here just call .save on the
        image. channel works as in cell_image.
        """
        from ..tiling.read_tiles import read_tiles

        mbtiles, channel, channels, note = self._background(channel, mbtiles)
        img_w, img_h, planes = _mbtiles_meta(mbtiles)
        if plane is None:
            plane, why = planes[len(planes) // 2], 'the middle plane of the stack'
        else:
            plane, why = int(plane), 'asked for'
            if plane not in planes:
                raise ValueError('no plane %s in %s, it has planes %s to %s'
                                 % (plane, mbtiles.name, planes[0], planes[-1]))
        if bbox is not None:
            x0, y0, x1, y1 = map(float, bbox)
            bbox = (max(x0, 0.0), max(y0, 0.0), min(x1, img_w), min(y1, img_h))
        im, scale = read_tiles(str(mbtiles), plane=plane, bbox=bbox, width=width)
        box = bbox or (0.0, 0.0, float(img_w), float(img_h))
        info = {
            'plane': plane,
            'plane_is': why,
            'planes_in_the_stack': [planes[0], planes[-1]],
            'image_size': [img_w, img_h],
            'bbox': [round(float(v), 1) for v in box],
            'scale': round(float(scale), 3),
            'to_place_a_point': 'a point (x, y) of the image is at ((x - bbox[0]) * scale, '
                                '(y - bbox[1]) * scale) in this picture',
            'background': channel,
            'backgrounds_in_this_run': channels,
            **({'background_note': note} if note else {}),
            'mbtiles': str(mbtiles),
            'image_is': 'stitched from the jpeg tile pyramid, a close visual copy of the '
                        'image, not the raw pixels. Fine to look at, not to measure',
        }
        return im, info

    def _has_containment(self):
        """Whether this run carries the inside_cell column, added Sept 2026."""
        shard = self.path.parent.parent / 'arrow_spots'
        files = sorted(shard.glob('*.feather')) if shard.is_dir() else []
        if not files:
            return False
        import pyarrow.feather as feather
        return 'inside_cell' in feather.read_table(files[0]).schema.names

    def _nb_contributions(self, c):
        """Per gene, per class negative binomial log-likelihood for one cell.

        The same arithmetic as compute_gene_loglikelihood_matrix, for one row:
        the expected count is the class mean scaled by the cell's area, by the gene
        efficiency and by the cell's theta, with SpotReg added on.
        """
        eta = np.asarray(self._meta('eta_bar', parse=True), dtype=np.float32)
        spot_reg = float(self._meta('SpotReg'))
        r_spot = float(self._meta('rSpot'))

        expected = np.einsum('gk,g,k->gk', c['scaled_means'], eta, c['theta_bar']) + spot_reg
        p = expected / (r_spot + expected)
        # the counts the class update scored, when the run saved them; gene_count is
        # after the spot update that follows it (see cell_to_cellType)
        counts = c.get('scored_gene_count', c['gene_count'])
        return counts[:, None] * np.log(p) + r_spot * np.log(1 - p)


# below this many soft counts a gene is taken as absent from the cell. Such a gene
# can still be strong evidence: a class that expresses it is penalised for the
# cell not holding any. The story must say so, or a reader (or an agent) sees a
# gene with zero counts listed as evidence and takes it for a mistake.
_ABSENT = 0.1


def shape_of(cov, mcr=None):
    """A covariance matrix as the shape of a cell, in numbers and in words.

    cov is d x d, d = 2 or 3, in pixels of the xy plane. mcr is the mean cell
    radius, the radius of the prior sphere, or None when the run does not say.
    """
    cov = np.asarray(cov, dtype=np.float64)
    d = cov.shape[0]
    axes = ['x', 'y', 'z'][:d]
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    radii = np.sqrt(np.maximum(vals[order], 0))
    vecs = vecs[:, order]
    # a zero radius cannot come out of the update, V0 is positive, but guard it
    # anyway: an infinity would not survive the trip through JSON
    elongation = float(radii[0] / radii[-1]) if radii[-1] > 0 else None

    if elongation is None:
        shape = 'degenerate, no spread along one axis'
    elif elongation < 1.1:
        shape = 'round'
    elif elongation < 1.5:
        shape = 'slightly elongated'
    else:
        shape = 'elongated'

    # the direction of the longest axis, up to sign, so make its biggest
    # component positive. Meaningless for a round cell
    long_axis = None
    long_axis_is = 'not given, the cell is round so no axis stands out'
    if elongation is not None and elongation >= 1.1:
        v = vecs[:, 0]
        i = int(np.argmax(np.abs(v)))
        v = v * np.sign(v[i])
        long_axis = {a: float(round(v[j], 3)) for j, a in enumerate(axes)}
        if abs(v[i]) >= 0.8:
            words = 'mostly along %s' % axes[i]
        else:
            # the runner up component; not argsort[-2], which on a tie gives i again
            j = [int(k) for k in np.argsort(-np.abs(v), kind='stable') if k != i][0]
            words = 'tilted, between %s and %s' % (axes[i], axes[j])
        long_axis_is = ('a unit vector, the direction of the longest axis, %s. z counts '
                        'in pixels of the xy plane, like Z in cellData' % words)

    out = {
        'shape': shape,
        'radii': [float(round(r, 3)) for r in radii],
        'radii_are': 'the standard deviation of the cell\'s Gaussian along each of its own '
                     'axes, longest first, in pixels of the xy plane. The 3 sigma '
                     'ellipsoid, cellData\'s sphere_scale, is three times these',
        'elongation': None if elongation is None else float(round(elongation, 3)),
        'elongation_is': 'longest radius over shortest; 1 is a sphere',
        'long_axis': long_axis,
        'long_axis_is': long_axis_is,
        'cov': [[float(round(x, 4)) for x in r] for r in cov],
        'cov_is': 'the covariance of the Gaussian, %s order, in pixels of the xy plane '
                  'squared' % ' '.join(axes),
    }
    if mcr is not None:
        out['radii_vs_mean_cell_radius'] = [float(round(r / mcr, 3)) for r in radii]
        out['radii_vs_mean_cell_radius_is'] = (
            'each radius over the mean cell radius, %.3f. The prior shape is a sphere '
            'of that radius, so 1, 1, 1 is a cell the data did not move away from it'
            % mcr)
        if np.allclose(radii, mcr, rtol=1e-3):
            out['shape'] = 'round, the prior sphere of the mean cell radius'
            out['note'] = ('the covariance is the prior one. Either the run was made '
                           'without the cell geometry update, or cell_cov_prior was '
                           'large enough for the prior to dominate')
    return out


def _present_absent(rows, n=3):
    """The top genes split into those the cell holds and those it lacks.

    Names for the present ones, whole rows for the absent ones, since the story
    quotes what a cell of the other class typically holds.
    """
    top = rows[:n]
    return ([g['gene'] for g in top if g['counts'] >= _ABSENT],
            [g for g in top if g['counts'] < _ABSENT])


def _about(x):
    """One decimal, rounded first so the JS port gets the same digits."""
    return '%.1f' % (round(x * 10) / 10)


def narrate_cell(e):
    """The story behind an explain_cell result, in plain words.

    Built from the numbers rather than left to the agent, so the story is the same
    whoever asks and cannot get the sign of a log-likelihood the wrong way round.
    Written for a reader who has never heard of a nat: the score differences are
    turned into odds ('about sixty to one') or into words ('overwhelmingly'), and the
    genes are ranked rather than numbered. The numbers themselves stay in the dict
    for anyone who wants them.
    """
    a, o = e['assigned'], e['compared_with']
    s = e['score']
    d = {
        'genes': s['gene_loglik']['assigned'] - s['gene_loglik']['compared'],
        'prior': s['log_prior']['assigned'] - s['log_prior']['compared'],
        'neighbours': s['spatial']['assigned'] - s['spatial']['compared'],
    }
    margin = sum(d.values())
    decider = max(d, key=lambda k: abs(d[k]))
    for_a = [g['gene'] for g in e['genes_favouring_assigned']]
    for_o = [g['gene'] for g in e['genes_favouring_compared']]
    a_present, a_absent = _present_absent(e['genes_favouring_assigned'])
    o_present, o_absent = _present_absent(e['genes_favouring_compared'])

    def p(x):
        return 'less than 0.01' if x < 0.005 else '%.2f' % x

    def names(gs, n=3):
        gs = gs[:n]
        return gs[0] if len(gs) == 1 else ', '.join(gs[:-1]) + ' and ' + gs[-1]

    # the absent genes and, in the same order, what a cell of the class that does
    # express them typically holds: 'Pcp4' and '3.6', 'Pcp4 and Car2' and '3.6 and 1.2'
    def absent(rows, key):
        return names([g['gene'] for g in rows]), names([_about(g[key]) for g in rows])

    out = []
    out.append('Cell %d was called %s, with probability %s. The closest alternative was %s, '
               'at %s.' % (e['cell'], a, p(e['prob_assigned']), o, p(e['prob_compared'])))
    out.append('pciSeq decides a cell\'s class from three things: how well its gene counts '
               'match what each class expresses according to the cell type definitions '
               '(the gene log-likelihood), how '
               'common each class is to begin with (the prior), and what the neighbouring '
               'cells were called (the spatial term). The class that comes out best '
               'overall wins.')

    # the genes, ranked, no numbers
    if d['genes'] > 0:
        if a_present:
            # 'typically' used to stand here; say where it comes from, the cell type
            # definitions, so it is not read as the run's own average (that is mean_in_*)
            out.append('The genes point to %s, %s. The strongest evidence comes from %s: the '
                       'cell holds these in amounts that fit what the cell type definitions '
                       'give for a %s cell, and not for a %s cell.'
                       % (a, _strength(d['genes']), names(a_present), a, o))
            if a_absent:
                gs, ms = absent(a_absent, 'mean_in_compared')
                out.append('%s %s the same way by %s absence: the cell holds almost none, '
                           'where cells this run called %s hold about %s on average.'
                           % (gs, 'counts' if len(a_absent) == 1 else 'count',
                              'its' if len(a_absent) == 1 else 'their', o, ms))
        else:
            gs, ms = absent(a_absent, 'mean_in_compared')
            out.append('The genes point to %s, %s. The strongest evidence is absence: the cell '
                       'holds almost no %s, where cells this run called %s hold about %s on '
                       'average.'
                       % (a, _strength(d['genes']), gs, o, ms))
        if o_present:
            out.append('A few genes, %s, look more like %s, but they are outweighed.'
                       % (names(o_present), o))
        if o_absent:
            gs, ms = absent(o_absent, 'mean_in_assigned')
            out.append('The near absence of %s also leans towards %s, since cells this run '
                       'called %s hold about %s on average, but not by enough.' % (gs, o, a, ms))
    else:
        if o_present:
            out.append('On its genes alone the cell looks more like %s, %s, mostly because of '
                       '%s.' % (o, _strength(-d['genes']), names(o_present)))
            if o_absent:
                gs, ms = absent(o_absent, 'mean_in_assigned')
                out.append('The near absence of %s says the same: cells this run called %s '
                           'hold about %s on average.' % (gs, a, ms))
        else:
            gs, ms = absent(o_absent, 'mean_in_assigned')
            out.append('On its genes alone the cell looks more like %s, %s, mostly through '
                       'absence: it holds almost no %s, where cells this run called %s hold '
                       'about %s on average.' % (o, _strength(-d['genes']), gs, a, ms))
        if a_present:
            out.append('The genes arguing for %s are %s.' % (a, names(a_present)))
        if a_absent:
            gs, ms = absent(a_absent, 'mean_in_compared')
            out.append('The near absence of %s argues for %s too, since cells this run '
                       'called %s hold about %s on average.' % (gs, a, o, ms))

    # the prior
    if abs(d['prior']) < 0.05:
        out.append('The prior treats the two classes alike.')
    else:
        who = a if d['prior'] > 0 else o
        out.append('The prior favours %s, because that class is more common to begin with.'
                   % who)

    # the neighbourhood
    if abs(d['neighbours']) < 0.4:
        out.append('The neighbouring cells make no real difference either way.')
    elif d['neighbours'] > 0:
        # the spatial term says which of the two classes the neighbourhood favours,
        # not how many neighbours are of it: espio cell 2413 got 'the neighbouring
        # cells are mostly Ndnf' with 2 of its 9 neighbours Ndnf
        out.append('The neighbourhood %s favours %s over %s, which %s the call.'
                   % (_mostly(d['neighbours']), a, o,
                      'strengthens' if d['genes'] > 0 else 'is what carries'))
    else:
        out.append('The neighbourhood favours %s, which counts against the call.' % o)

    # what settled it
    if decider == 'genes':
        out.append('So the genes settled it%s.'
                   % (', and the neighbourhood agreed' if d['neighbours'] > 0.4 else ''))
    elif decider == 'neighbours' and d['genes'] <= 0:
        out.append('So this call is the neighbourhood overruling the genes: on its genes '
                   'alone the cell would have been called %s, but surrounded by %s cells '
                   'the balance comes out for %s, %s.' % (o, a, a, _strength(margin)))
    elif decider == 'neighbours':
        out.append('So the neighbourhood settled it. The genes agreed, but only mildly; the '
                   'surrounding cells made the difference.')
    else:
        out.append('So the prior settled it.')
    return ' '.join(out)


def narrate_spot(e):
    """The story behind an explain_spot result, in plain words.

    Same idea as narrate_cell: computed from the numbers so it is the same whoever
    asks, no units in the output, differences said as odds or as a word. The shape
    follows the spot page of the docs: where the spot is against how well its gene
    fits each cell, and which of the two carried the call.
    """
    gene, spot = e['gene'], e['spot']
    cells = [c for c in e['candidates'] if c['cell'] != 'background']
    bg = next(c for c in e['candidates'] if c['cell'] == 'background')
    by_prob = sorted(cells, key=lambda c: -c['prob'])
    winner = by_prob[0]
    to_bg = e['assigned_to'] == 'background'

    def expr(c):
        return c['class expression'] + c['cell scale'] + c['cell-gene scale']

    def p(x):
        return 'less than 0.01' if x < 0.005 else '%.2f' % x

    def sure(x):
        return ('confidently' if x > 0.9 else 'fairly confidently' if x > 0.6
                else 'narrowly' if x > 0.4 else 'with no clear winner')

    out = []
    # the call
    if to_bg:
        out.append('Spot %d is a %s spot. It was assigned to the background, with probability '
                   '%s, meaning the model takes it for a misread rather than a read from any '
                   'cell. The nearest cell, %d, gets %s.'
                   % (spot, gene, p(bg['prob']), winner['cell'], p(winner['prob'])))
    else:
        others = ', '.join('cell %d (%s)' % (c['cell'], p(c['prob'])) for c in by_prob[1:3])
        out.append('Spot %d is a %s spot. It was assigned to cell %d, %s, with probability %s. '
                   'The next candidates are %s, and the chance it is a misread is %s.'
                   % (spot, gene, winner['cell'], sure(winner['prob']), p(winner['prob']),
                      others, p(bg['prob'])))

    # how a spot is scored
    out.append('pciSeq weighs each nearby cell on two things: how close the spot is to the '
               'cell\'s centroid (the Gaussian fit), and how well a %s spot fits that cell, '
               'which combines whether the cell\'s class expresses %s (class expression), '
               'whether the cell holds more reads overall than its class predicts (cell '
               'scale), and whether it already holds more %s than its class predicts '
               '(cell-gene scale). The background is scored on how often %s spots turn out '
               'to be misreads. The best total wins.' % (gene, gene, gene, gene))

    if to_bg:
        out.append('Here no cell scores well enough: cell %d is the nearest, and its class '
                   'is %s, but %s does not fit it well, and %s misreads are common enough in '
                   'this run for the background to win.'
                   % (winner['cell'], winner['class'], gene, gene))
        return ' '.join(out)

    # position
    by_dist = sorted(cells, key=lambda c: -c['Gaussian fit'])
    nearest = by_dist[0]
    if nearest['cell'] == winner['cell']:
        rivals = [c for c in by_dist[1:3]]
        out.append('Cell %d is the nearest candidate: %s on position alone.'
                   % (winner['cell'], ' and '.join(
                       '%s over cell %d' % (_strength(winner['Gaussian fit'] - c['Gaussian fit']), c['cell'])
                       for c in rivals)))
    else:
        out.append('Cell %d is not the nearest: cell %d is closer, %s on position alone.'
                   % (winner['cell'], nearest['cell'],
                      _strength(nearest['Gaussian fit'] - winner['Gaussian fit'])))

    # expression
    runner = by_prob[1] if len(by_prob) > 1 else None
    classes = {c['class'] for c in cells}
    if len(classes) == 1:
        # a shared most likely class is not enough: the term averages over the class
        # probabilities, so cells of the same class can still differ on it. Name the
        # measured best rather than guess the cause
        ce = [c['class expression'] for c in cells]
        if max(ce) - min(ce) <= 0.2:
            out.append('Every candidate is a %s cell, and on class alone %s fits them about '
                       'equally; what separates them is how much each already holds.'
                       % (winner['class'], gene))
        else:
            best = max(cells, key=lambda c: c['class expression'])
            out.append('Every candidate is a %s cell, but their class probabilities differ, so '
                       'on class alone %s fits cell %d best.' % (winner['class'], gene, best['cell']))
    elif runner and runner['class'] != winner['class']:
        # compared between the two, never read on its own: the term is
        # log(class mean x Inefficiency + SpotReg), below zero whenever the class mean
        # is under about 9 counts, so its sign does not say whether a class expresses
        # the gene (spot 18377 on espio: both classes express Cnr1, both terms < 0)
        d_ce = winner['class expression'] - runner['class expression']
        if abs(d_ce) <= 0.2:
            out.append('Cell %d is a %s cell and cell %d a %s cell, and the two classes '
                       'express %s about equally.'
                       % (winner['cell'], winner['class'], runner['cell'], runner['class'], gene))
        else:
            out.append('Cell %d is a %s cell, whose class expresses %s %s than the class of '
                       'cell %d, %s.'
                       % (winner['cell'], winner['class'], gene, 'more' if d_ce > 0 else 'less',
                          runner['cell'], runner['class']))
    elif runner:
        out.append('Cell %d is a %s cell, and so is cell %d.'
                   % (winner['cell'], winner['class'], runner['cell']))
    else:
        out.append('Cell %d is a %s cell.' % (winner['cell'], winner['class']))
    if runner:
        pieces = []
        d_sc = winner['cell scale'] - runner['cell scale']
        d_cg = winner['cell-gene scale'] - runner['cell-gene scale']
        if abs(d_sc) > 0.2:
            pieces.append('cell %d holds more reads overall than its class predicts'
                          % (winner['cell'] if d_sc > 0 else runner['cell']))
        if abs(d_cg) > 0.2:
            pieces.append('cell %d already holds more %s than its class predicts'
                          % (winner['cell'] if d_cg > 0 else runner['cell'], gene))
        if pieces:
            out.append('Between the top two, %s.' % ' and '.join(pieces))

    # the inside cell bonus, when it is switched on
    bonus = [c for c in cells if c.get('bonus', 0)]
    if bonus:
        out.append('The spot\'s pixel lies inside the mask of cell %d, which adds a bonus for '
                   'it.' % bonus[0]['cell'])

    # verdict, winner against the runner up
    if runner:
        d_pos = winner['Gaussian fit'] - runner['Gaussian fit']
        d_exp = expr(winner) - expr(runner)
        if d_pos > 0.4 and d_exp > 0.4:
            out.append('So cell %d is both the nearer and the better fit for the gene, and '
                       'the call is clear.' % winner['cell'])
        elif d_pos > 0.4:
            out.append('So cell %d is the better fit for the gene, %s, but cell %d is closer, '
                       '%s, and distance carries the call.'
                       % (runner['cell'], _strength(-d_exp) if d_exp < -0.4 else 'slightly',
                          winner['cell'], _strength(d_pos)))
        elif d_exp > 0.4:
            out.append('So cell %d is closer, %s, but cell %d fits the gene better, %s, and '
                       'expression carries the call.'
                       % (runner['cell'], _strength(-d_pos) if d_pos < -0.4 else 'slightly',
                          winner['cell'], _strength(d_exp)))
        else:
            out.append('So the two are close on both counts, and the call is a narrow one.')
    if bg['prob'] > 0.1:
        out.append('There is a real chance, %s, that the spot is a misread.' % p(bg['prob']))
    return ' '.join(out)


def _strength(d):
    """A log-likelihood difference as something a reader can picture: the odds, e**d
    to one, in round numbers while the number still means anything, and a word once
    it does not. Meant to close a sentence: 'the genes point to X, about 60 to one'."""
    import math
    if d < 0.4:
        return 'only just'
    # a d above 10 is past 20,000 to one (e**10 is about 22,000), where the words
    # no longer change, so it is answered before taking exp: a cell against a far
    # off class (2413 vs Oligo on espio) has a d so big that math.exp raises
    # OverflowError
    if d > 10:
        return 'overwhelmingly, beyond any doubt'
    odds = math.exp(d)
    if odds < 3:
        return 'slightly'
    if odds < 20:
        return 'about %d to one' % round(odds)
    if odds < 100:
        return 'about %d to one' % (round(odds / 5) * 5)
    if odds < 20000:
        return 'about %s to one' % format(int(round(odds, -2)), ',')
    return 'overwhelmingly, beyond any doubt'


def _mostly(d):
    """How strongly the neighbourhood favours one class over the other, from the size of
    the spatial term difference. Not how many neighbours are of it, the term does not
    say that."""
    return 'strongly' if d > 8 else 'clearly' if d > 2 else 'slightly'


def _tsv(v):
    """A float the way the tsv files print it: three decimals.

    Only used on the rebuild path. The viewer files were already written from the
    rounded frame, so this is a no-op on them in practice, it is here so the rebuilt
    row cannot come out with more digits than the real file would have had.
    """
    return round(float(v), 3)


def _channel_name(mbtiles):
    """The name a background image goes by: the 'name' in its metadata, as the viewer
    labels its radio buttons, or the file name when there is none."""
    from ..tiling.read_tiles import _metadata
    con = sqlite3.connect('file:%s?mode=ro' % mbtiles, uri=True)
    try:
        name = (_metadata(con).get('name') or '').strip()
    finally:
        con.close()
    return name or Path(mbtiles).stem


def _mbtiles_meta(mbtiles):
    """Width and height of the original image and the planes, off the mbtiles
    metadata. Older files may have only plane_count, or neither for a 2D image."""
    from ..tiling.read_tiles import _metadata
    con = sqlite3.connect('file:%s?mode=ro' % mbtiles, uri=True)
    try:
        meta = _metadata(con)
    finally:
        con.close()
    if meta.get('planes'):
        planes = sorted(int(p) for p in meta['planes'].split(','))
    else:
        planes = list(range(int(meta.get('plane_count', 1))))
    return int(meta['width']), int(meta['height']), planes


# ------------------------------------------------------------- arithmetic

# calculate() walks the parsed expression itself and only knows these. Never eval:
# an expression is text from the model, and eval would run anything.
_CALC_FUNCS = {'exp': math.exp, 'log': math.log, 'log10': math.log10,
               'sqrt': math.sqrt, 'abs': abs, 'round': round}
_CALC_BINOPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
                ast.Div: operator.truediv, ast.Pow: operator.pow}
_CALC_MAX_LEN = 500       # characters, a sum of a few dozen numbers fits easily
_CALC_MAX_POWER = 1000    # |exponent|, so 10**10**10 cannot run the machine out


def calculate(expression):
    """Arithmetic for the agent, so it does not have to do it in its head.

    Numbers, + - * / ** and brackets, and exp, log (natural), log10, sqrt, abs and
    round. For example '15.29 + 9.8 + 3.1' to add up some gene differences, or
    'exp(3.2)' to turn a log-likelihood difference into odds. Anything else, a
    name, an attribute, a keyword argument, is refused with a message saying so.
    """
    text = str(expression)
    if len(text) > _CALC_MAX_LEN:
        raise ValueError('expression longer than %d characters' % _CALC_MAX_LEN)
    try:
        tree = ast.parse(text, mode='eval')
    except SyntaxError:
        raise ValueError('cannot read %r as arithmetic' % text)
    try:
        result = _calc(tree.body)
    except (OverflowError, ZeroDivisionError) as e:
        raise ValueError('%s: %s' % (text, e))
    except ValueError as e:
        # math domain errors, log(0) and the like, come here too
        raise ValueError('%s: %s' % (text, e))
    return {'expression': text, 'result': result}


def _calc(node):
    if isinstance(node, ast.Constant) and type(node.value) in (int, float):
        return node.value
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        v = _calc(node.operand)
        return -v if isinstance(node.op, ast.USub) else v
    if isinstance(node, ast.BinOp) and type(node.op) in _CALC_BINOPS:
        a, b = _calc(node.left), _calc(node.right)
        if isinstance(node.op, ast.Pow):
            if abs(b) > _CALC_MAX_POWER:
                raise ValueError('exponent larger than %d' % _CALC_MAX_POWER)
            # in floats, so an integer power cannot grow a number with a million digits
            return math.pow(a, b)
        return _CALC_BINOPS[type(node.op)](a, b)
    if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id in _CALC_FUNCS and not node.keywords):
        return _CALC_FUNCS[node.func.id](*[_calc(x) for x in node.args])
    raise ValueError('only numbers, + - * / ** and %s are allowed, not %s'
                     % (', '.join(_CALC_FUNCS), ast.dump(node)[:60]))


def _polygon_area(xs, ys):
    """Shoelace area of a closed polygon, in pixels."""
    x, y = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def _draw_outlines(im, scale, box, rows, colour, width, fill_alpha):
    """Draw outline rows (x_list, y_list in image pixels) onto a read_tiles window,
    filled see-through so the nuclei still show."""
    from PIL import Image, ImageDraw
    over = Image.new('RGBA', im.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(over)
    for r in rows:
        pts = [((x - box[0]) * scale, (y - box[1]) * scale)
               for x, y in zip(r['x_list'], r['y_list'])]
        if len(pts) < 3:
            continue
        d.polygon(pts, fill=colour + (int(255 * fill_alpha),))
        d.line(pts + [pts[0]], fill=colour + (255,), width=width, joint='curve')
    im.paste(Image.alpha_composite(im.convert('RGBA'), over).convert('RGB'), (0, 0))


def _find_db(path):
    """Accept the run folder, anything above it, or the database itself."""
    path = Path(path).expanduser()
    if path.is_file():
        return path
    hits = sorted(path.rglob('diagnostics.db'))
    if not hits:
        raise FileNotFoundError('no diagnostics.db under %s' % path)
    if len(hits) > 1:
        raise ValueError('%d diagnostics.db files under %s, point at one of them:\n  %s'
                         % (len(hits), path, '\n  '.join(str(h) for h in hits[:10])))
    return hits[0]


def open_run(path):
    """Open a finished run. The entry point every other tool goes through."""
    return Run(path)
