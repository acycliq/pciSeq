"""Tools for asking about a run while it is still running.

The finished-run tools (tools.py) read diagnostics.db, which does not exist until
the run ends. These read the VarBayes object in memory, so they answer the
questions only a live run can raise: how far has it got, is it converging, what is
this cell right now, what changed since the last iteration.

They are plain functions over a VarBayes, as tools.py is plain functions over a
Run: the live viewer's chat is a thin wrapper, and a notebook can call them just as
well while a fit runs in another thread.

Three rules, the same as the finished-run tools:
  1. cell labels in and out are the labels of the segmentation, never internal.
  2. counts are soft, weighted by assignment probability, and every answer says so.
  3. row 0 of the arrays is the background pseudocell, never reported.

Everything here is a snapshot: the loop keeps iterating while the agent reads, so
two answers in one conversation can come from different iterations. Every answer
carries the iteration it was taken at, and the agent is told to quote it.
"""
import numpy as np

from ..core.utils.cell_utils import fetch_label, reversed_label_map, to_internal

# counts below this are dropped from per gene listings, as in tools.py
_COUNT_TOL = 0.001


class Live:
    """A run in progress, opened for questions.

    Parameters
    ----------
    varbayes : VarBayes
        The model being fitted. Held by reference, never copied: the point is to
        see the current state, not a stale one.

    Notes
    -----
    The history of what each cell was called is not kept by the model, so
    `remember()` has to be called once per iteration for the movement tools to
    have anything to compare against. The live viewer calls it from the same
    callback that streams the class updates.
    """

    def __init__(self, varbayes):
        self.vb = varbayes
        self._prev_assigned = None      # argmax per cell at the previous iteration
        self._prev_iter = None
        self._flips = None              # how many times each cell has changed class
        self._reverse = None            # internal label -> segmentation label

    # ------------------------------------------------------------ plumbing

    @property
    def _cfg(self):
        return self.vb.config or {}

    def _label_map(self):
        return self._cfg.get('label_map')

    def to_internal(self, label):
        try:
            label = int(label)
        except (TypeError, ValueError):
            raise ValueError('%r is not a cell label' % (label,))
        try:
            row = to_internal(label, self._label_map())
        except KeyError:
            raise KeyError('no cell %s in this segmentation' % label)
        if row == 0:
            raise ValueError('cell 0 is the background pseudocell, not a cell')
        # without a label_map the label IS the row, so nothing has checked it yet.
        # Before the first update there is nothing to check against; the caller's
        # own guard then says the run has not got that far.
        arr = self.vb.cells.classProb
        if arr is None:
            arr = self.vb.cells.geneCount
        if arr is not None and not 0 < row < len(arr):
            raise KeyError('no cell %s in this run, it has %d' % (label, len(arr) - 1))
        return row

    def to_external(self, row):
        # to_external flips the whole map on every call, which is fine once but not
        # per row; keep the reversed map, the label_map cannot change mid run
        lm = self._label_map()
        if not lm:
            return int(row)
        if self._reverse is None:
            self._reverse = reversed_label_map(lm)
        return int(fetch_label(int(row), self._reverse))

    def _class_prob(self):
        p = self.vb.cells.classProb
        if p is None:
            raise RuntimeError('the run has not reached its first class update yet, '
                               'so there is nothing to report about the cells')
        return p

    @property
    def class_names(self):
        return np.asarray(self.vb.cells.class_names)

    @property
    def gene_panel(self):
        return np.asarray(self.vb.genes.gene_panel)

    def remember(self):
        """Note what every cell is called right now, so the next call to
        `changed_cells` can say what moved. Called once per iteration."""
        p = self.vb.cells.classProb
        if p is None:
            return
        now = np.argmax(p, axis=1)
        if self._flips is None:
            self._flips = np.zeros(len(now), dtype=np.int64)
        elif self._prev_assigned is not None and len(self._prev_assigned) == len(now):
            self._flips += (now != self._prev_assigned)
        self._prev_assigned = now.copy()
        self._prev_iter = self.vb.iter_num

    # --------------------------------------------------------------- tools

    def progress(self):
        """Where the run has got to, and whether it is settling."""
        delta = list(self.vb.iter_delta or [])
        cfg = self._cfg
        tol = cfg.get('CellCallTolerance')
        out = {
            'iteration': self.vb.iter_num,
            'max_iter': cfg.get('max_iter'),
            'converged': bool(self.vb.has_converged),
            'tolerance': tol,
            'delta': delta[-1] if delta else None,
            'delta_is': 'the largest change in a spot assignment probability over the '
                        'last iteration; the loop stops when it falls below the tolerance',
            'delta_history': [float(d) for d in delta[-20:]],
            'history_is': 'the last 20 iterations, oldest first',
        }
        if len(delta) >= 4:
            recent = np.asarray(delta[-4:], dtype=float)
            falling = bool(recent[-1] < recent[0])
            out['trend'] = ('falling' if falling else 'not falling') + \
                           ' over the last 4 iterations'
            # a delta that stops falling while well above the tolerance is the
            # plateau Dimitris knows well; say it plainly rather than guess
            if not falling and tol is not None and recent[-1] > tol:
                out['note'] = ('the change is not falling and is still above the '
                               'tolerance, so the loop is not settling at the moment')
        return out

    def cell(self, label):
        """What a cell is called at this iteration, and on what evidence."""
        row = self.to_internal(label)
        p = self._class_prob()[row]
        counts = self.vb.cells.geneCount[row]
        order = np.argsort(-p, kind='stable')
        top = np.argsort(-counts, kind='stable')[:10]
        out = {
            'cell': int(label),
            'iteration': self.vb.iter_num,
            'total_counts': float(counts.sum()),
            'counts_are': 'soft, weighted by the spot assignment probabilities, as they '
                          'stand at this iteration',
            'classes': [{'class': str(self.class_names[k]), 'prob': float(p[k])}
                        for k in order[:5]],
            'top_genes': [{'gene': str(self.gene_panel[g]), 'counts': float(counts[g])}
                          for g in top if counts[g] > _COUNT_TOL],
        }
        if self._prev_assigned is not None and row < len(self._prev_assigned):
            was = int(self._prev_assigned[row])
            now = int(order[0])
            out['changed_since_last_iteration'] = was != now
            if was != now:
                out['was'] = str(self.class_names[was])
            if self._flips is not None:
                out['times_it_has_changed_class'] = int(self._flips[row])
        return out

    def class_counts(self):
        """How many cells each class holds at this iteration, hard and soft."""
        p = self._class_prob()[1:]          # row 0 is the background
        hard = np.bincount(np.argmax(p, axis=1), minlength=p.shape[1])
        soft = p.sum(axis=0)
        rows = [{'class': str(name), 'cells': int(hard[k]), 'soft': float(soft[k])}
                for k, name in enumerate(self.class_names)]
        zero = [r for r in rows if r['class'] == 'Zero']
        rest = sorted([r for r in rows if r['class'] != 'Zero'], key=lambda r: -r['cells'])
        return {
            'iteration': self.vb.iter_num,
            'n_cells': int(len(p)),
            'classes': zero + rest,
            'cells_is': 'hard: the number of cells whose most probable class it is now',
            'soft_is': 'the class probability summed over the cells',
        }

    def changed_cells(self, n=20):
        """The cells that changed class in the last iteration, and the ones that
        keep changing."""
        if self._prev_assigned is None:
            raise RuntimeError('nothing has been remembered yet, so there is no previous '
                               'iteration to compare against')
        now = np.argmax(self._class_prob(), axis=1)
        if len(now) != len(self._prev_assigned):
            raise RuntimeError('the number of cells changed since the last snapshot')
        moved = np.where(now != self._prev_assigned)[0]
        moved = moved[moved > 0]
        rows = [{'cell': self.to_external(i),
                 'from': str(self.class_names[self._prev_assigned[i]]),
                 'to': str(self.class_names[now[i]]),
                 'prob': float(self._class_prob()[i, now[i]])}
                for i in moved[:n]]
        out = {
            'iteration': self.vb.iter_num,
            'compared_with_iteration': self._prev_iter,
            'n_changed': int(len(moved)),
            'shown': len(rows),
            'cells': rows,
            'cells_are': 'the cells whose most probable class differs from the previous '
                         'snapshot, at most n of them',
        }
        if self._flips is not None:
            busy = np.argsort(-self._flips, kind='stable')
            busy = [i for i in busy if i > 0 and self._flips[i] > 1][:n]
            out['keeps_changing'] = [
                {'cell': self.to_external(i), 'times': int(self._flips[i]),
                 'now': str(self.class_names[now[i]])} for i in busy]
            out['keeps_changing_is'] = ('the cells that have changed class most often '
                                        'since the chat started watching; a cell flipping '
                                        'every iteration is not settling')
        return out

    def scale_factors(self):
        """The factors the model has learned so far: eta per gene, theta per cell."""
        genes = self.vb.genes
        eta = np.asarray(genes.eta_bar) if getattr(genes, 'eta_bar', None) is not None else None
        theta = self.vb.cells.theta_bar
        out = {'iteration': self.vb.iter_num}
        if eta is not None:
            order = np.argsort(-eta, kind='stable')
            out['eta'] = {
                'is': 'eta_bar, the posterior mean of eta, the gene\'s detection efficiency '
                      'relative to the Inefficiency setting; the reference '
                      'expression of every class is multiplied by it',
                'mean': float(eta.mean()),
                'highest': [{'gene': str(self.gene_panel[g]), 'eta': float(eta[g])}
                            for g in order[:5]],
                'lowest': [{'gene': str(self.gene_panel[g]), 'eta': float(eta[g])}
                           for g in order[-5:]],
            }
        if theta is not None:
            p = self._class_prob()
            scalar = np.einsum('ck,ck->c', p, theta)[1:]
            out['theta'] = {
                'is': "the cell scale factor, averaged over the class probabilities; it "
                      "scales a class's expected counts to the cell's total",
                'mean': float(scalar.mean()),
                'min': float(scalar.min()),
                'max': float(scalar.max()),
            }
        return out
