"""Fill in the answer key of espio_demo.yaml from the run, through the agent's tools.

No number is written in the yaml. Each case names its cell, spot or gene, and this
asks the same tools the agent uses for the facts a correct answer must hold, then
prints every question with its key filled in. A case picked for a property, a near
tie for instance, is flagged when the run no longer has it, so it can be pointed
at another cell. No language model is called, it costs nothing.

    PYTHONPATH=. python tests/eval/answer_key.py      (from the repo, not the installed copy)
"""
import os
import sys

import pyarrow.feather as feather
import yaml

from pciSeq.src.core.utils.geometry import plane_of
from pciSeq.src.mcp.tools import open_run

HERE = os.path.dirname(os.path.abspath(__file__))
r2 = lambda x: '%.2f' % x


def cell_confident(r, s):
    e, top = r.explain_cell(s['cell']), r.cell(s['cell'])['classes']
    return ({'cell': s['cell'], 'assigned': e['assigned'], 'prob': r2(e['prob_assigned']),
             'runner_up': top[1]['class']},
            e['prob_assigned'] > 0.9, 'the cell is no longer a confident call')


def cell_close_call(r, s):
    top = r.cell(s['cell'])['classes']
    return ({'cell': s['cell'], 'assigned': top[0]['class'], 'prob': r2(top[0]['prob']),
             'runner_up': top[1]['class'], 'runner_up_prob': r2(top[1]['prob'])},
            top[0]['prob'] - top[1]['prob'] < 0.1, 'the top two classes are no longer close')


def cell_why_not(r, s):
    e = r.explain_cell(s['cell'], s['compared'])
    return ({'cell': s['cell'], 'compared': s['compared']},
            e['prob_compared'] < 0.01, 'the compared class is no longer far off')


def spot_split(r, s):
    e = r.explain_spot(s['spot'])
    c = sorted((x for x in e['candidates'] if x['cell'] != 'background'), key=lambda x: -x['prob'])
    return ({'spot': s['spot'], 'gene': e['gene'], 'cell': c[0]['cell'], 'prob': r2(c[0]['prob']),
             'other': c[1]['cell'], 'other_prob': r2(c[1]['prob'])},
            c[0]['prob'] - c[1]['prob'] < 0.1, 'the spot is no longer split between two cells')


def rtheta(r, s):
    return {'rTheta': r.run_info()['settings']['rTheta']}, True, ''


def converged(r, s):
    i = r.run_info()
    return {'converged': 'yes' if i['converged'] else 'no', 'iterations': i['iterations']}, True, ''


def zero_class(r, s):
    z = [c for c in r.class_counts()['classes'] if c['class'] == 'Zero'][0]
    return {'cells': z['cells']}, True, ''


def tac1_counts(r, s):
    g = r.gene(s['gene'])
    return {'spots': g['total_spots'], 'misreads': g['spots_called_misread']}, True, ''


def plane_floor(r, s):
    path = os.path.dirname(os.path.dirname(r.path))
    cells = feather.read_table(os.path.join(path, 'arrow_cells', 'cells_shard_000.feather'),
                               columns=['cell_id', 'Z']).to_pandas().set_index('cell_id')
    z, vs = float(cells.loc[s['cell'], 'Z']), r.run_info()['voxel_size']
    exact = z * vs[0] / vs[2]
    return ({'cell': s['cell'], 'plane': plane_of(z, vs), 'rounded': round(exact),
             'plane_exact': r2(exact), 'z': r2(z)},
            round(exact) != plane_of(z, vs), 'floor and round now agree, the case no longer tells them apart')


def tac1_efficiency(r, s):
    eta = r.gene(s['gene'])['eta']
    return {'eta': r2(eta), 'inverted': r2(1 / eta)}, True, ''


def no_such_cell(r, s):
    try:
        r.cell(s['cell'])
        return {'cell': s['cell']}, False, 'the run has this cell after all'
    except KeyError:
        return {'cell': s['cell']}, True, ''


KEYS = {f.__name__: f for f in (cell_confident, cell_close_call, cell_why_not, spot_split, rtheta,
                                converged, zero_class, tac1_counts, plane_floor, tac1_efficiency,
                                no_such_cell)}


def fill(text, values):
    return text.format(**values) if isinstance(text, str) else text


def main():
    spec = yaml.safe_load(open(os.path.join(HERE, 'espio_demo.yaml')))
    r = open_run(os.path.expanduser(spec['run']))
    stale = []
    for case in spec['cases']:
        values, holds, why = KEYS[case['id']](r, case.get('subject', {})) if case['id'] in KEYS else ({}, True, '')
        print('\n[%s]%s %s' % (case['id'], ' (viewer only)' if case.get('viewer_only') else '',
                                fill(case['question'], values)))
        for m in case.get('must', []):
            print('   must:     ' + fill(str(m), values))
        for m in case.get('must_not', []):
            print('   must not: ' + fill(str(m), values))
        if not holds:
            stale.append(case['id'])
            print('   !! ' + why)
    if stale:
        sys.exit('\n%d cases need another subject: %s' % (len(stale), ', '.join(stale)))


if __name__ == '__main__':
    main()
