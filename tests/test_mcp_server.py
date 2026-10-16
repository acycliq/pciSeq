"""The MCP wrapper round the tools.

tools.py is tested on its own in test_mcp_tools.py. This checks the thin layer on
top: every tool is registered, calls go through, and the errors tools.py raises on
purpose reach the client with their message intact. That last one matters most,
mcp hides the message of any exception it does not recognise, so without the
wrapping the agent would see 'Error executing tool' and nothing else.

Skipped when mcp is not installed, it is an optional extra.
"""
import asyncio
import json

import numpy as np
import pytest

mcp = pytest.importorskip('mcp')

from pciSeq.src.mcp import server as srv
from tests.test_label_identifiers import _run

TOOLS = {'open_run', 'cell', 'cell_shape', 'explain_cell', 'explain_spot', 'cell_counts', 'spots_in_cell',
         'spots_of_cell', 'cell_row', 'spot_row', 'docs', 'run_info', 'cell_image',
         'plane_image', 'gene', 'theta', 'gamma', 'spot', 'neighbours', 'class_counts',
         'class_prior', 'find_cells', 'metadata', 'calculate', 'list_source', 'read_source'}


def test_calculate_answers_and_refusals_reach_the_agent():
    """calculate needs no run, and when it refuses the agent must see why."""
    assert call('calculate', expression='15.29 + 9.8')['result'] == pytest.approx(25.09)
    msg = error_of('calculate', expression='__import__("os")')
    assert 'only numbers' in msg


def call(name, **args):
    """Call a tool the way a client does, and hand back the parsed JSON."""
    res = asyncio.run(srv.server.call_tool(name, args))
    assert not res.is_error
    return json.loads(res.content[0].text)


def error_of(name, **args):
    """The message an agent would see when a tool refuses."""
    from mcp.server.mcpserver.exceptions import ToolError
    with pytest.raises(ToolError) as exc:
        asyncio.run(srv.server.call_tool(name, args))
    return str(exc.value)


@pytest.fixture(scope='module')
def run_folder(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('mcp_server')
    _run(np.random.default_rng(3), tmp, save=True, max_iter=40, tol=0.001)
    return tmp


def test_every_tool_is_registered():
    names = {t.name for t in asyncio.run(srv.server.list_tools())}
    assert names == TOOLS


def test_every_tool_has_a_description_the_agent_can_use():
    for t in asyncio.run(srv.server.list_tools()):
        assert t.description and len(t.description) > 40, t.name


def test_calls_go_through(run_folder):
    s = call('open_run', path=str(run_folder))
    assert s['cells'] == 16

    e = call('explain_cell', label=105)
    assert e['assigned'] != e['compared_with']

    c = call('cell_counts', label=105)
    assert 'soft' in c['counts_are']

    h = call('spots_in_cell', label=105)
    assert 'no probabilities' in h['spots_are']

    a = call('spots_of_cell', label=105)
    b = call('spots_of_cell', label=105, min_prob=0.0001)
    assert a['n_spots'] <= b['n_spots']


def test_the_image_tools_come_back_as_pictures(run_folder, tmp_path):
    """The agent has to get an actual png it can look at, plus the facts as text."""
    import base64
    from tests.test_mcp_images import _grey_mbtiles
    mb = tmp_path / 'grey.mbtiles'
    _grey_mbtiles(mb)
    call('open_run', path=str(run_folder))
    out = tmp_path / 'cell105.png'
    res = asyncio.run(srv.server.call_tool(
        'cell_image', {'label': 105, 'mbtiles': str(mb), 'width': 300, 'save_as': str(out)}))
    assert not res.is_error
    assert [c.type for c in res.content] == ['image', 'text']
    png = base64.b64decode(res.content[0].data)
    assert res.content[0].mime_type == 'image/png' and png[:4] == b'\x89PNG'
    info = json.loads(res.content[1].text)
    assert info['cell'] == 105 and info['saved_as'] == str(out)
    assert out.read_bytes() == png

    res = asyncio.run(srv.server.call_tool('plane_image', {'mbtiles': str(mb), 'width': 100}))
    assert [c.type for c in res.content] == ['image', 'text']
    assert json.loads(res.content[1].text)['bbox'] == [0.0, 0.0, 80.0, 80.0]


def test_the_agent_is_told_to_ask_which_background(tmp_path):
    """Two backgrounds and no channel: the agent must get the names and the
    instruction to ask, not a bare 'Error executing tool'."""
    from tests.test_mcp_images import _grey_mbtiles
    _run(np.random.default_rng(3), tmp_path, save=True)
    call('open_run', path=str(tmp_path))
    vd = next(tmp_path.rglob('diagnostics.db')).parent.parent
    _grey_mbtiles(vd / 'dapi.mbtiles', name='DAPI')
    _grey_mbtiles(vd / 'gcamp.mbtiles', name='GCaMP')
    msg = error_of('plane_image')
    assert 'DAPI' in msg and 'GCaMP' in msg and 'Ask the user' in msg
    res = asyncio.run(srv.server.call_tool('plane_image', {'channel': 'GCaMP', 'width': 50}))
    assert json.loads(res.content[1].text)['background'] == 'GCaMP'


def test_the_message_survives_when_a_tool_refuses(run_folder):
    """The whole reason server.py wraps the tools."""
    call('open_run', path=str(run_folder))
    assert 'no cell 999999' in error_of('cell', label=999999)
    assert 'background' in error_of('cell_counts', label=0)
    assert 'NoSuchGene' in error_of('cell_counts', label=105, gene='NoSuchGene')


def test_nothing_works_before_open_run():
    srv._run = None
    assert 'open_run' in error_of('cell', label=105)


def test_open_run_on_nothing_says_so(tmp_path):
    assert 'no diagnostics.db' in error_of('open_run', path=str(tmp_path))


# ------------------------------------------------------- packaging

def test_the_module_is_shipped_and_launchable():
    """Two ways a plain pip install can leave the server unreachable: the package
    missing from find_packages, or no console script so the only way in is the
    internal module path."""
    import ast
    import pathlib
    from setuptools import find_packages

    repo = pathlib.Path(__file__).resolve().parents[1]
    assert 'pciSeq.src.mcp' in find_packages(str(repo))

    setup_py = (repo / 'setup.py').read_text()
    assert 'pciseq-mcp = pciSeq.src.mcp.server:main' in setup_py
    assert callable(srv.main)


def test_mcp_is_an_extra_not_a_core_dependency():
    """A plain install must not drag in a web server, a jwt stack and a pydantic
    floor for the sake of a feature most users never touch. Same call spikelab
    makes."""
    import pathlib
    import runpy
    import unittest.mock as mock

    repo = pathlib.Path(__file__).resolve().parents[1]
    with mock.patch('setuptools.setup') as setup:
        runpy.run_path(str(repo / 'setup.py'), run_name='not_main')
    kwargs = setup.call_args.kwargs

    core = ' '.join(kwargs['install_requires'])
    assert 'mcp' not in core, 'mcp belongs in extras_require, not install_requires'

    extras = kwargs['extras_require']
    assert extras['mcp'] == ['mcp>=2']
    # no [all] while mcp is the only extra: it would mean the same thing under a
    # second name, which only confuses. Bring it back with the second extra.
    assert 'all' not in extras


def test_missing_mcp_still_says_something_useful():
    """mcp ships with pciSeq now, so this only fires on an old install, but the
    message should still say what to do."""
    import pathlib
    src = (pathlib.Path(__file__).resolve().parents[1]
           / 'pciSeq/src/mcp/server.py').read_text()
    assert 'pciSeq_3d[mcp]' in src


def test_the_source_tools_read_github_at_the_commit_that_made_the_run(run_folder, monkeypatch):
    from pciSeq.src.mcp import source
    call('open_run', path=str(run_folder))
    commit = call('run_info')['commit']
    asked = []

    def fake_github(url):
        asked.append(url)
        if url == source.RAW + commit + '/pciSeq/src/core/main.py':
            return 'line one\nline two\nline three\n'
        raise FileNotFoundError(url)
    monkeypatch.setattr(source, 'fetch_text', fake_github)

    out = call('read_source', path='core/main.py', start_line=2)
    assert out['path'] == 'pciSeq/src/core/main.py' and out['start_line'] == 2
    assert out['commit'] == commit and 'made this run' in out['source_is']
    assert asked == [source.RAW + commit + '/pciSeq/src/core/main.py']
    # the refusals reach the agent with their message
    assert 'not a path inside the repo' in error_of('read_source', path='../setup.py')
    assert 'only the python source' in error_of('read_source', path='pciSeq/src/tiling/README.md')
    assert 'never pushed' in error_of('read_source', path='core/nope.py')


def test_the_core_persona_is_short_and_copied_word_for_word():
    """The server sends CORE, not the full INSTRUCTIONS: Claude Code cuts a server's
    instructions at 4,096 characters. Every sentence of CORE must be in the full
    persona word for word, so a rewording there that is not made here too fails, and
    CORE must stay under the cut with room to spare."""
    import re
    from pciSeq.src.mcp.server_instructions import CORE, INSTRUCTIONS
    from pciSeq.src.mcp.server import server
    norm = lambda t: re.sub(r'\s+', ' ', t).strip()
    full = norm(INSTRUCTIONS)
    sentences = [x for p in CORE.split('\n\n') for x in re.split(r'(?<=[.;])\s+(?=[A-Z"])', p)]
    assert len(sentences) > 20
    missing = [x for x in sentences if norm(x) not in full]
    assert not missing, missing
    assert len(CORE) < 3900, len(CORE)
    # plain ASCII, so there is no doubt how the 4,096 is counted
    assert CORE.isascii()
    assert server.instructions == CORE


def test_no_tool_description_is_long_enough_to_be_cut():
    """Claude Code cuts a tool description at 4,096 characters when it loads it up
    front, 16,384 through tool search. Anything past the cut never reaches the model,
    so every description stays well under the smaller one."""
    for t in asyncio.run(srv.server.list_tools()):
        assert len(t.description) < 4000, (t.name, len(t.description))
        assert t.description.isascii(), t.name


def test_the_tool_rules_are_the_personas_own_words():
    """Each tool in TOOL_RULES carries passages cut out of INSTRUCTIONS, for an agent
    that gets only CORE. If the persona is reworded and a start or an end no longer
    matches, tool_rules raises; and every passage must really sit on its tool."""
    import re
    from pciSeq.src.mcp.server_instructions import TOOL_RULES, tool_rules, INSTRUCTIONS
    full = re.sub(r'\s+', ' ', INSTRUCTIONS)
    listed = {t.name: re.sub(r'\s+', ' ', t.description) for t in asyncio.run(srv.server.list_tools())}
    for name, passages in TOOL_RULES.items():
        for start, end in passages:
            i = full.find(start)
            assert i >= 0 and full.find(end, i) >= 0, (name, start, end)
        assert tool_rules(name) in listed[name], name

