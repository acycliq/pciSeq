// Ask the test questions through the viewer's own chat and save the answers.
//
// The questions come from answer_key.py, filled in from the run:
//
//   PYTHONPATH=. python tests/eval/answer_key.py --json questions.json
//   node tests/eval/run_chat_panel.js questions.json answers.json
//
// It starts pciSeq Viewer on the run named in espio_demo.yaml, asks each question
// as a new conversation with the chat's own runTurn (the same code the chat panel
// runs, same persona, same tools), saves the answer and the tools it called, and
// closes the viewer. The answers are then marked against the must / must not lines
// that come with each question.
//
// COSTS API CREDIT: it uses the provider, model and key saved in the viewer's chat
// settings, one conversation per question. Never part of pytest or npm test.
//
// Needs node 22 or later, and pciSeq Viewer at VIEWER below (or the --viewer
// option). Close any open viewer first: this one uses the same settings folder and
// the same MCP port. The viewer saves the folder it opens in its config, so the
// config is copied before and put back after.
//
//   --viewer DIR   the pciSeq_viewer repo (default ~/dev/javascript/pciSeq_viewer)
//   --data DIR     the viewer_data folder (default the run of espio_demo.yaml)

const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawn } = require('child_process');

const args = process.argv.slice(2);
const opt = (name, fallback) => {
  const i = args.indexOf(name);
  return i >= 0 ? args.splice(i, 2)[1] : fallback;
};
const VIEWER = opt('--viewer', path.join(os.homedir(), 'dev/javascript/pciSeq_viewer'));
const DATA = opt('--data', path.join(os.homedir(), 'pciseq_runs/espio_demo/pciSeq/data/viewer_data'));
const [questionsFile, answersFile] = args;
if (!questionsFile || !answersFile) {
  console.error('usage: node run_chat_panel.js questions.json answers.json [--viewer DIR] [--data DIR]');
  process.exit(1);
}
const MAIN_PORT = 9230;     // the node inspector of the viewer's main process
const PAGE_PORT = 9334;     // the devtools port of its window
const CHAT = path.join(VIEWER, 'electron', 'chat.js');
const CONFIG = path.join(os.homedir(), '.config', 'pciSeq_viewer', 'config.json');

const sleep = ms => new Promise(r => setTimeout(r, ms));

// run an expression over a devtools websocket and hand back its value
async function evaluate(port, expr, pick = list => list[0]) {
  const list = await (await fetch(`http://127.0.0.1:${port}/json`)).json();
  const ws = new WebSocket(pick(list).webSocketDebuggerUrl);
  return new Promise((resolve, reject) => {
    ws.onopen = () => ws.send(JSON.stringify({ id: 1, method: 'Runtime.evaluate',
      params: { expression: expr, awaitPromise: true, returnByValue: true, timeout: 600000 } }));
    ws.onmessage = ev => {
      const r = JSON.parse(ev.data);
      if (r.id !== 1) return;
      ws.close();
      if (r.result.exceptionDetails) reject(new Error(JSON.stringify(r.result.exceptionDetails).slice(0, 600)));
      else resolve(r.result.result.value);
    };
    ws.onerror = reject;
  });
}
const inMain = expr => evaluate(MAIN_PORT, expr);
const inPage = expr => evaluate(PAGE_PORT, expr, list => list.find(t => t.type === 'page'));

async function waitFor(what, test, seconds) {
  for (let i = 0; i < seconds; i++) {
    try { if (await test()) return; } catch (e) { /* not up yet */ }
    await sleep(1000);
  }
  throw new Error(`the viewer did not get ${what} within ${seconds} s`);
}

(async () => {
  const cases = JSON.parse(fs.readFileSync(questionsFile, 'utf8'));
  const configCopy = fs.existsSync(CONFIG) ? fs.readFileSync(CONFIG) : null;
  const viewer = spawn(path.join(VIEWER, 'node_modules', '.bin', 'electron'),
    ['.', `--inspect=${MAIN_PORT}`, `--remote-debugging-port=${PAGE_PORT}`, '--remote-allow-origins=*', `--data=${DATA}`],
    { cwd: VIEWER, stdio: 'ignore' });
  try {
    await waitFor('its inspector open', async () => (await fetch(`http://127.0.0.1:${MAIN_PORT}/json`)).ok, 60);
    // the spots stream in after the window opens; questions about genes need them
    await waitFor('the run loaded', async () => (await inPage('String(window.appState?.geneDataMap?.size || 0)')) !== '0', 300);
    const model = JSON.parse(await inMain(
      `JSON.stringify((({provider, model}) => ({provider, model}))(process.mainModule.require(${JSON.stringify(CHAT)}).getSettings()))`));
    console.log('model', model.provider, model.model, '|', cases.length, 'questions');

    const results = [];
    for (const c of cases) {
      const t0 = Date.now();
      const row = { ...c };
      try {
        const got = JSON.parse(await inMain(
          `process.mainModule.require(${JSON.stringify(CHAT)})` +
          `.runTurn([{ role: 'user', content: ${JSON.stringify(c.question)} }])` +
          `.then(r => JSON.stringify({ text: r.text, tools: r.messages.flatMap(m => Array.isArray(m.content)` +
          ` ? m.content.filter(b => b.type === 'tool_use').map(b => ({ name: b.name, input: b.input })) : []) }))`));
        Object.assign(row, { answer: got.text, called: got.tools });
      } catch (e) {
        row.error = e.message;
      }
      row.seconds = Math.round((Date.now() - t0) / 1000);
      results.push(row);
      console.log(c.id, row.error ? 'ERROR ' + row.error.slice(0, 150)
                                  : 'tools: ' + row.called.map(t => t.name).join(', '), row.seconds + ' s');
      // written after every question, so a run cut short keeps what it has
      fs.writeFileSync(answersFile, JSON.stringify({ channel: 'chat_panel', model, date: new Date().toISOString(),
                                                     data: DATA, results }, null, 1));
    }
  } finally {
    viewer.kill();
    if (configCopy) fs.writeFileSync(CONFIG, configCopy);
  }
})().catch(e => { console.error(e.message); process.exit(1); });
