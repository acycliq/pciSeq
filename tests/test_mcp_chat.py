"""The live viewer's chat: one turn of conversation.

The model is stubbed, so these check what we control: that the tool loop runs the
live tools and hands the results back in the API's shape, and that each service
gets its key the way it wants it. The settings file itself is test_chat_settings.

"""
import json

import pytest

from pciSeq.src.mcp import chat
from pciSeq.src.mcp.live import Live


@pytest.fixture
def cfg(tmp_path):
    return tmp_path / 'pciseq' / 'chat.json'


@pytest.fixture
def live(minimal_varbayes):
    vb = minimal_varbayes
    vb.initialise_state()
    vb.iter_num = 0
    vb.geneCount_upd()
    vb.rho_upd()
    vb.eta_upd()
    vb.theta_upd()
    vb.gamma_upd()
    vb.cell_to_cellType()
    vb.iter_delta.append(0.3)
    lv = Live(vb)
    lv.remember()
    return lv


# -------------------------------------------------------------- one turn

class FakeResponse:
    def __init__(self, payload, status=200):
        self._payload = payload
        self.status_code = status
        self.text = json.dumps(payload)

    def json(self):
        return self._payload


def _reply(content, stop='end_turn'):
    return {'id': 'm', 'type': 'message', 'role': 'assistant', 'model': 'x',
            'content': content, 'stop_reason': stop,
            'usage': {'input_tokens': 1, 'output_tokens': 1}}


def test_a_turn_runs_the_tool_and_hands_the_result_back(live):
    sent = []
    calls = []

    def post(url, **kw):
        sent.append((url, kw))
        if len(sent) == 1:
            return FakeResponse(_reply(
                [{'type': 'tool_use', 'id': 't1', 'name': 'progress', 'input': {}}], 'tool_use'))
        return FakeResponse(_reply([{'type': 'text', 'text': 'it is on iteration 0'}]))

    out = chat.run_turn(live, [{'role': 'user', 'content': 'how far along is it?'}],
                        on_event=calls.append,
                        settings={'provider': 'anthropic', 'api_key': 'sk-test',
                                  'model': 'claude-sonnet-5', 'base_url': None,
                                  'key_from_env': False},
                        post=post)
    assert out['text'] == 'it is on iteration 0'
    url, kw = sent[0]
    assert url == 'https://api.anthropic.com/v1/messages'
    assert kw['headers']['x-api-key'] == 'sk-test'
    assert 'authorization' not in kw['headers']
    assert [t['name'] for t in kw['json']['tools']] == [t['name'] for t in chat.TOOLS]
    assert 'live viewer' in kw['json']['system'] and 'neuroscientist' in kw['json']['system']
    # the tool ran against the live model and its answer went back as a tool_result
    result = sent[1][1]['json']['messages'][2]['content'][0]
    assert result['type'] == 'tool_result' and result['tool_use_id'] == 't1'
    assert json.loads(result['content'])['iteration'] == 0
    assert [e['type'] for e in calls] == ['tool_call', 'tool_result', 'text', 'done']


def test_a_compatible_service_gets_a_bearer_token(live):
    def post(url, **kw):
        assert url == 'https://api.z.ai/api/anthropic/v1/messages'
        assert kw['headers']['authorization'] == 'Bearer zai-key'
        assert 'x-api-key' not in kw['headers']
        return FakeResponse(_reply([{'type': 'text', 'text': 'ok'}]))

    chat.run_turn(live, [{'role': 'user', 'content': 'hi'}],
                  settings={'provider': 'zai', 'api_key': 'zai-key', 'model': 'glm-5.2',
                            'base_url': 'https://api.z.ai/api/anthropic', 'key_from_env': False},
                  post=post)


def test_a_tool_that_refuses_reaches_the_model(live):
    sent = []

    def post(url, **kw):
        sent.append(kw)
        if len(sent) == 1:
            return FakeResponse(_reply(
                [{'type': 'tool_use', 'id': 't1', 'name': 'cell', 'input': {'label': 99999}}],
                'tool_use'))
        return FakeResponse(_reply([{'type': 'text', 'text': 'no such cell'}]))

    chat.run_turn(live, [{'role': 'user', 'content': 'tell me about cell 99999'}],
                  settings={'provider': 'anthropic', 'api_key': 'sk', 'model': 'm',
                            'base_url': None, 'key_from_env': False},
                  post=post)
    answer = json.loads(sent[1]['json']['messages'][2]['content'][0]['content'])
    assert 'no cell 99999' in answer['error']


def test_the_model_s_own_error_message_is_passed_on(live):
    def post(url, **kw):
        return FakeResponse({'error': {'message': 'credit balance is too low'}}, status=400)

    with pytest.raises(RuntimeError, match='credit balance is too low'):
        chat.run_turn(live, [{'role': 'user', 'content': 'hi'}],
                      settings={'provider': 'anthropic', 'api_key': 'sk', 'model': 'm',
                                'base_url': None, 'key_from_env': False},
                      post=post)


def test_no_key_is_refused_before_anything_is_sent(live, cfg):
    with pytest.raises(RuntimeError, match='no API key'):
        chat.run_turn(live, [{'role': 'user', 'content': 'hi'}],
                      settings={'provider': 'anthropic', 'api_key': '', 'model': 'm',
                                'base_url': None, 'key_from_env': False},
                      post=lambda *a, **k: pytest.fail('nothing should be sent'))


def test_the_persona_is_the_one_from_the_python_server():
    from pciSeq.src.mcp.server_instructions import INSTRUCTIONS
    prompt = chat.system_prompt()
    assert prompt.startswith(INSTRUCTIONS)
    assert 'reads the model AS IT IS NOW' in prompt
    assert 'no diagnostics.db yet' in prompt
