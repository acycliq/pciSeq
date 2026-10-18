"""The live viewer's chat: one turn of conversation about the running fit.

The page sends a question, this runs the tool loop against the live model and
sends back the answer. The API key stays here, in the server process, and never
reaches the browser (chat_settings.py).

The loop is the standard one: send the messages with the tool list; if the model
answers with tool_use blocks, run each against the Live object, hand the results
back as tool_result blocks, and go round again; stop when it answers with text
only. Every step is reported through `on_event` so the page can show what was
looked up rather than only the final answer.

Any service speaking Anthropic's Messages API works, so the model can be Claude or
GLM, as in the desktop viewer. One POST per round, with requests, which pciSeq
already depends on through flask; no SDK.
"""
import json
import logging

import requests

from . import chat_settings
from .live import Live
from .server_instructions import INSTRUCTIONS

logger = logging.getLogger(__name__)

MAX_TOKENS = 2048
MAX_TOOL_ROUNDS = 8
TIMEOUT = 120

# What is only true of a run in progress. The shared persona (INSTRUCTIONS) says
# how to explain a pciSeq result; this says what is different when the answer is a
# moving target.
LIVE_SYSTEM = '\n'.join([
    'You are inside the pciSeq live viewer, the window that opens while a run is',
    'still being fitted. The user is watching the cells change colour as the loop',
    'iterates. Every tool here reads the model AS IT IS NOW, so two answers in one',
    'conversation can come from different iterations: each answer carries the',
    'iteration it was taken at, and you should say which iteration you are talking',
    'about whenever it matters.',
    '',
    'Nothing here is final. A cell that is Zero at iteration 12 may not be at',
    'iteration 40, so speak about what the model thinks so far rather than what the',
    'run found, and never call a result final while the loop is running. progress',
    'says how far it has got and whether the change per iteration is still falling;',
    'changed_cells says what moved since the last iteration and which cells keep',
    'flipping, which is the question behind "is this settling".',
    '',
    'There is no diagnostics.db yet, so the tools that read a finished run are not',
    'here: no spot level scores, no pictures, no cellData. When the user asks for',
    'something only a finished run can answer, say that it will be available in the',
    'pciSeq viewer once the run has saved its output.',
])

# the tools the model is offered, built from the Live methods
TOOLS = [
    {
        'name': 'progress',
        'description':
            'Where the run has got to: the iteration, the largest change in a spot '
            'assignment over the last iteration (delta) with its recent history, the '
            'tolerance the loop stops at, and whether the change is still falling. Use '
            'it for "how far along is it", "is it converging", "how many iterations '
            'left", "is it settling".',
        'input_schema': {'type': 'object', 'properties': {}},
    },
    {
        'name': 'cell',
        'description':
            'What a cell is called at this iteration and on what evidence: its class '
            'probabilities, its top genes by soft count, whether it changed class since '
            'the last iteration and how many times it has changed since the chat '
            'started watching. label is the cell number of the segmentation.',
        'input_schema': {
            'type': 'object',
            'properties': {'label': {'type': 'integer',
                                     'description': 'The cell label, as in the segmentation.'}},
            'required': ['label'],
        },
    },
    {
        'name': 'class_counts',
        'description':
            'How many cells each class holds at this iteration, hard (the cells whose '
            'most probable class it is) and soft (the class probability summed over the '
            'cells). Zero first, the rest by size. Use it for "how many cells are Zero '
            'now", "which classes are growing".',
        'input_schema': {'type': 'object', 'properties': {}},
    },
    {
        'name': 'changed_cells',
        'description':
            'What moved since the last iteration: the cells whose most probable class '
            'changed, from which class to which, and separately the cells that keep '
            'changing. A cell flipping every iteration is not settling. Use it for '
            '"what changed", "is anything still moving", "which cells are unstable".',
        'input_schema': {
            'type': 'object',
            'properties': {'n': {'type': 'integer', 'description': 'How many to list, default 20.'}},
        },
    },
    {
        'name': 'scale_factors',
        'description':
            'The factors the model has learned so far: eta per gene (the gene '
            'efficiency, with the highest and lowest few) and theta per cell (the cell '
            'scale factor, as a range). Use it for "what is eta now", "have the scale '
            'factors settled".',
        'input_schema': {'type': 'object', 'properties': {}},
    },
    {
        'name': 'docs',
        'description':
            'The pciSeq documentation. query returns the paragraphs holding those '
            'words, best first, each with its page and heading; page returns one whole '
            'page. The search matches words literally, so when your words may not be '
            'the ones the docs use, read a page instead. Use this before answering any '
            'question about how pciSeq works, a term, or a setting such as rTheta, '
            'mrf_beta or Inefficiency.',
        'input_schema': {
            'type': 'object',
            'properties': {
                'query': {'type': 'string', 'description': 'A few words, for example "rTheta".'},
                'page': {'type': 'string', 'description': 'A page, for example "the-model/scale-factors.md".'},
                'n': {'type': 'integer', 'description': 'How many paragraphs, for a query; default 5.'},
            },
        },
    },
    {
        'name': 'settings',
        'description':
            'The settings this run was started with (mrf_beta, rTheta, Inefficiency, '
            'nNeighbors, CellCallTolerance, max_iter and the rest), and the size of the '
            'run. Use it for "what settings is this using", "what is the tolerance".',
        'input_schema': {'type': 'object', 'properties': {}},
    },
    {
        'name': 'list_source',
        'description':
            'List a folder of the pciSeq source code, the code that is running this '
            'fit. Use it to find the file a question is about; the model code is under '
            'pciSeq/src/core.',
        'input_schema': {
            'type': 'object',
            'properties': {'dir': {'type': 'string', 'description': 'Folder path, "" for the top.'}},
        },
    },
    {
        'name': 'read_source',
        'description':
            'Read a file of the pciSeq source code, the code that is running this fit, '
            'with line numbers. Reach for it only when a question needs the actual '
            'code, after the docs; never claim to run it. Long files come back in '
            'slices of up to 400 lines, give start_line to read on. When you cite it, '
            'give the file and the line. A path from the code map page works as it is.',
        'input_schema': {
            'type': 'object',
            'properties': {
                'path': {'type': 'string', 'description': 'For example "pciSeq/src/core/main.py".'},
                'start_line': {'type': 'integer', 'description': 'First line to return, default 1.'},
            },
            'required': ['path'],
        },
    },
]


def system_prompt():
    return INSTRUCTIONS + '\n\n' + LIVE_SYSTEM


def call_tool(live, name, args):
    """Run one tool against the live model. Errors come back as a dict the model
    can read, rather than as an exception, so it can tell the user what happened."""
    try:
        if name == 'progress':
            return live.progress()
        if name == 'cell':
            return live.cell(args['label'])
        if name == 'class_counts':
            return live.class_counts()
        if name == 'changed_cells':
            return live.changed_cells(n=args.get('n', 20))
        if name == 'scale_factors':
            return live.scale_factors()
        if name == 'settings':
            cfg = dict(live.vb.config or {})
            cfg.pop('label_map', None)             # 25k entries, not for the model
            return {
                'settings': cfg,
                'cells': int(len(live.vb.cells.geneCount) - 1),
                'genes': int(len(live.gene_panel)),
                'classes': int(len(live.class_names)),
            }
        if name == 'docs':
            from . import docs as docs_mod
            # a run being fitted has no saved docs yet, so these are this machine's,
            # the pages of the pciSeq doing the fit
            page = args.get('page')
            if page:
                try:
                    text = docs_mod.read_page(page)
                except KeyError:
                    return {'error': 'no page %s' % page,
                            'pages': docs_mod.list_pages()}
                return {'page': page, 'title': docs_mod.page_title(text), 'text': text}
            hits = docs_mod.search_docs(args.get('query', ''), n=args.get('n', 5))
            return {'query': args.get('query', ''), 'hits': hits,
                    'docs_are': docs_mod.source()}
        if name == 'list_source':
            from . import source as source_mod
            return source_mod.list_source(args.get('dir', ''))
        if name == 'read_source':
            from . import source as source_mod
            return source_mod.read_source(args['path'], args.get('start_line', 1))
        return {'error': 'unknown tool %s' % name}
    except Exception as e:
        return {'error': '%s' % e}


def run_turn(varbayes, messages, on_event=None, settings=None, post=None):
    """One turn: the model, its tool calls, and the answer.

    varbayes : the model being fitted, held by reference
    messages : the conversation so far, in the API's shape, ending with the user's
               new message. Returned extended, so the next turn has the context.
    on_event : called with {'type': ...} as the turn proceeds: text, tool_call,
               tool_result, done, error. The page shows these.
    post     : for the tests, in place of requests.post
    """
    s = settings or chat_settings.load()
    if not s['api_key']:
        raise RuntimeError('no API key. Set one on the chat settings, or put it in %s.'
                           % (chat_settings.PROVIDERS[s['provider']]['env'] or 'the settings'))
    say = on_event or (lambda ev: None)
    live = Live(varbayes) if not isinstance(varbayes, Live) else varbayes
    url = (s['base_url'] or 'https://api.anthropic.com').rstrip('/') + '/v1/messages'
    # anthropic takes the key as x-api-key, the compatible services as a bearer
    headers = {'anthropic-version': '2023-06-01', 'content-type': 'application/json'}
    if s['provider'] == 'anthropic':
        headers['x-api-key'] = s['api_key']
    else:
        headers['authorization'] = 'Bearer ' + s['api_key']

    history = list(messages)
    final_text = ''
    send = post or requests.post

    for _ in range(MAX_TOOL_ROUNDS + 1):
        r = send(url, headers=headers, timeout=TIMEOUT, json={
            'model': s['model'],
            'max_tokens': MAX_TOKENS,
            'system': system_prompt(),
            'tools': TOOLS,
            'messages': history,
        })
        if r.status_code != 200:
            # the model's own message is the useful part, not the status code
            try:
                detail = r.json().get('error', {}).get('message', r.text[:300])
            except ValueError:
                detail = r.text[:300]
            raise RuntimeError('the model answered %d: %s' % (r.status_code, detail))
        res = r.json()
        content = res.get('content', [])
        history.append({'role': 'assistant', 'content': content})

        texts = [b['text'] for b in content if b.get('type') == 'text']
        if texts:
            final_text = '\n'.join(texts)
            say({'type': 'text', 'text': final_text})

        uses = [b for b in content if b.get('type') == 'tool_use']
        if res.get('stop_reason') != 'tool_use' or not uses:
            break

        results = []
        for u in uses:
            say({'type': 'tool_call', 'name': u['name'], 'input': u.get('input', {})})
            out = call_tool(live, u['name'], u.get('input', {}))
            say({'type': 'tool_result', 'name': u['name'], 'is_error': 'error' in out})
            results.append({'type': 'tool_result', 'tool_use_id': u['id'],
                            'content': json.dumps(out, default=float)})
        history.append({'role': 'user', 'content': results})

    say({'type': 'done'})
    return {'text': final_text, 'messages': history}
