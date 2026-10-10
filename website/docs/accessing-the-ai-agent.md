---
description: The four ways to access the AI agent that explains a run, and what each one uses.
# standard column and margins; only the diagram runs wide, as on the scaling factors demo
pageClass: demo-wide
---

<script setup>
import { ref, computed } from 'vue'

// the path of a question: which step is picked, 0 is all of them
const step = ref(0)
const steps = [
  { name: 'All', says: 'The path of one question. Dashed arrows cross the network to the language model; solid arrows stay on the user\'s machine.' },
  { name: '1 Question', says: 'The question is sent to the language model with the list of functions: the name, the description and the arguments of each.' },
  { name: '2 Function to call', says: 'The language model picks one from the list and replies with its name and arguments, here explain_cell for cell 5016.' },
  { name: '3 Function runs', says: 'The function is run on the user\'s machine and reads the saved files.' },
  { name: '4 Result', says: 'The result of the function is sent to the language model.' },
  { name: '5 Answer', says: 'The language model writes the answer from that result.' },
]
const lit = (...on) => step.value === 0 || on.includes(step.value)

// requirements: one selection per row, and the lower block lists what that
// selection needs. A row that does not apply to the selection fades.
const when = ref('after')        // 'running' | 'after'
const where = ref('assistant')   // 'live' | 'panel' | 'assistant'
const tool = ref('claude')       // 'claude' | 'gemini' | 'codex'
const open = ref(true)           // the run is open in pciSeq Viewer
const pickWhen = (v) => {
  when.value = v
  if (v === 'running') where.value = 'live'
  else if (where.value === 'live') where.value = 'assistant'
}
const pickWhere = (v) => {
  if (v === 'live' && when.value !== 'running') return
  if (v !== 'live' && when.value === 'running') return
  where.value = v
}
// the summary table jumps straight to one result
const go = (k) => {
  when.value = k === 'live' ? 'running' : 'after'
  where.value = k === 'live' ? 'live' : k === 'panel' ? 'panel' : 'assistant'
  if (k === 'viewer') open.value = true
  if (k === 'mcp') open.value = false
}
const result = computed(() =>
  where.value === 'live' ? 'live' : where.value === 'panel' ? 'panel' : open.value ? 'viewer' : 'mcp')
const toolName = { claude: 'Claude Code', gemini: 'Gemini CLI', codex: 'Codex' }
const register = {
  viewer: {
    claude: 'claude mcp add --transport http pciseq-viewer http://127.0.0.1:8317/mcp',
    gemini: 'gemini mcp add --transport http pciseq-viewer http://127.0.0.1:8317/mcp',
    codex:  'codex mcp add pciseq-viewer --url http://127.0.0.1:8317/mcp',
  },
  mcp: {
    claude: 'claude mcp add --scope user pciseq -- pciseq-mcp',
    gemini: 'gemini mcp add --scope user pciseq pciseq-mcp',
    codex:  'codex mcp add pciseq -- pciseq-mcp',
  },
}
const copied = ref('')
const copy = (text) => {
  navigator.clipboard.writeText(text).then(() => {
    copied.value = text
    setTimeout(() => { if (copied.value === text) copied.value = '' }, 1500)
  })
}
const plain = 'pip install git+https://github.com/acycliq/pciSeq_3d.git@dev_3d'
const extra = 'pip install "pciSeq_3d[mcp] @ git+https://github.com/acycliq/pciSeq_3d.git@dev_3d"'
</script>

# Accessing the AI agent

pciSeq saves with its results, for every cell and every spot, the terms from which its
assignment was computed. An agent, a language model equipped with functions that read
those files, explains an assignment from them: why a cell was assigned to its class, why
a spot was assigned to one cell rather than another, what value a setting took.

There are four ways to put a question to it. The diagram shows what each one uses.
Choosing a way keeps its parts and fades the rest; the person marks where the question
is typed.

<DemoFrame src="/agent_access.html" wide :min="420"
           title="The ways to ask about a run" />

| Way | When | Reads | Requires |
| :--- | :--- | :--- | :--- |
| An AI assistant with `pciseq-mcp` | once the run has completed | the saved files, in the results folder | see <a href="#requirements" @click="go('mcp')">Requirements</a> |
| An AI assistant connected to [pciSeq Viewer](./viewer/overview.md) | once the run has completed | the saved files, as loaded in pciSeq Viewer | see <a href="#requirements" @click="go('viewer')">Requirements</a> |
| The chat panel in pciSeq Viewer | once the run has completed | the saved files, as loaded in pciSeq Viewer | see <a href="#requirements" @click="go('panel')">Requirements</a> |
| The chat of the [live viewer](./api/live-viewer.md) | during the run | the fit in progress | see <a href="#requirements" @click="go('live')">Requirements</a> |

## Requirements

The requirements depend on the stage of the run, the interface and, for an AI assistant,
on the assistant and on whether the results are open in pciSeq Viewer. The lower block lists
them for the selection. An AI assistant is a program on the user's machine that sends
questions to a language model and supports the Model Context Protocol (MCP), such as
Claude Code, Gemini CLI or Codex.

<div class="pick" role="group" aria-label="Requirements">

<div class="pick-row">
  <div class="pick-label">Stage</div>
  <div class="pick-opts">
    <button type="button" :class="{ on: when === 'running' }" :aria-pressed="when === 'running'" @click="pickWhen('running')">during the run</button>
    <button type="button" :class="{ on: when === 'after' }" :aria-pressed="when === 'after'" @click="pickWhen('after')">after the run</button>
  </div>
</div>

<div class="pick-row">
  <div class="pick-label">Interface</div>
  <div class="pick-opts">
    <button type="button" :class="{ on: where === 'live' }" :disabled="when !== 'running'" :aria-pressed="where === 'live'" @click="pickWhere('live')">live viewer chat</button>
    <button type="button" :class="{ on: where === 'panel' }" :disabled="when === 'running'" :aria-pressed="where === 'panel'" @click="pickWhere('panel')">pciSeq Viewer chat panel</button>
    <button type="button" :class="{ on: where === 'assistant' }" :disabled="when === 'running'" :aria-pressed="where === 'assistant'" @click="pickWhere('assistant')">AI assistant</button>
  </div>
</div>

<div class="pick-row" :class="{ off: where !== 'assistant' }">
  <div class="pick-label">Assistant</div>
  <div class="pick-opts">
    <button v-for="(name, k) in toolName" :key="k" type="button" :class="{ on: tool === k }" :disabled="where !== 'assistant'" :aria-pressed="tool === k" @click="tool = k">{{ name }}</button>
  </div>
</div>

<div class="pick-row" :class="{ off: where !== 'assistant' }">
  <div class="pick-label">Use pciSeq Viewer</div>
  <div class="pick-opts">
    <button type="button" :class="{ on: open }" :disabled="where !== 'assistant'" :aria-pressed="open" @click="open = true">yes</button>
    <button type="button" :class="{ on: !open }" :disabled="where !== 'assistant'" :aria-pressed="!open" @click="open = false">no</button>
  </div>
</div>

<div class="pick-out">

<div class="pick-row">
  <div class="pick-label">Installation</div>
  <div class="pick-v">
    <template v-if="result === 'live'"><div><a href="./installation">pciSeq</a>, without the <code>mcp</code> extra. The live viewer and its chat are part of it.</div><div class="pick-cmd"><code>{{ plain }}</code><button type="button" class="pick-copy" :class="{ done: copied === plain }" @click="copy(plain)" aria-label="Copy"><template v-if="copied === plain"><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M20 6 9 17l-5-5"/></svg></template><template v-else><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg></template></button></div></template>
    <template v-else-if="result === 'mcp'"><div><a href="./installation">pciSeq</a>, with the <code>mcp</code> extra, which <code>pciseq-mcp</code> needs.</div><div class="pick-cmd"><code>{{ extra }}</code><button type="button" class="pick-copy" :class="{ done: copied === extra }" @click="copy(extra)" aria-label="Copy"><template v-if="copied === extra"><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M20 6 9 17l-5-5"/></svg></template><template v-else><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg></template></button></div></template>
    <template v-else><a href="./viewer/overview">pciSeq Viewer</a>, from its <a href="https://github.com/acycliq/pciSeq_viewer/releases/latest">releases page</a>.</template>
  </div>
</div>

<div class="pick-row" :class="{ off: result === 'live' || result === 'panel' }">
  <div class="pick-label">Registration</div>
  <div class="pick-v">
    <template v-if="result === 'live' || result === 'panel'"></template>
    <template v-else-if="result === 'viewer'">{{ toolName[tool] }} is registered with pciSeq Viewer, once, by adding the address of pciSeq Viewer to {{ toolName[tool] }}:<div class="pick-cmd"><code>{{ register.viewer[tool] }}</code><button type="button" class="pick-copy" :class="{ done: copied === register.viewer[tool] }" @click="copy(register.viewer[tool])" aria-label="Copy"><template v-if="copied === register.viewer[tool]"><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M20 6 9 17l-5-5"/></svg></template><template v-else><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg></template></button></div></template>
    <template v-else>The <code>pciseq-mcp</code> program is registered with {{ toolName[tool] }}, once. {{ toolName[tool] }} starts it when needed.<div class="pick-cmd"><code>{{ register.mcp[tool] }}</code><button type="button" class="pick-copy" :class="{ done: copied === register.mcp[tool] }" @click="copy(register.mcp[tool])" aria-label="Copy"><template v-if="copied === register.mcp[tool]"><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M20 6 9 17l-5-5"/></svg></template><template v-else><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg></template></button></div></template>
  </div>
</div>

<div class="pick-row" :class="{ off: result !== 'mcp' }">
  <div class="pick-label">Note</div>
  <div class="pick-v">
    <template v-if="result === 'mcp'">The first question of a {{ toolName[tool] }} session gives the results folder; later questions in the same session do not.<div class="pick-cmd"><code>Open &lt;pciSeq results folder&gt; and explain why cell 100 is assigned to its class.</code><button type="button" class="pick-copy" :class="{ done: copied === 'Open <pciSeq results folder> and explain why cell 100 is assigned to its class.' }" @click="copy('Open <pciSeq results folder> and explain why cell 100 is assigned to its class.')" aria-label="Copy"><template v-if="copied === 'Open <pciSeq results folder> and explain why cell 100 is assigned to its class.'"><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M20 6 9 17l-5-5"/></svg></template><template v-else><svg viewBox="0 0 24 24" width="16" height="16" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg></template></button></div></template>
  </div>
</div>

<div class="pick-row">
  <div class="pick-label">API key</div>
  <div class="pick-v">
    <template v-if="result === 'live'">Required.</template>
    <template v-else-if="result === 'panel'">Required.</template>
    <template v-else>Not required.</template>
  </div>
</div>

</div>
</div>

An AI assistant without pciSeq Viewer covers a machine where it is not installed or not
used, for example a remote machine reached over ssh. The two ways to connect an
assistant are independent: each has its own command and neither requires the other.
Any other application that supports MCP connects the same way; only the form of the
registration command differs, see the [MCP server](./api/mcp-server.md) page.

## The path of a question

In every way a question is answered in the same five steps. The language model does
not read the saved files. It names a function, and the function is run on the
user's machine. Choosing a step keeps its parts and fades the rest.

The language model has no knowledge of these functions of its own. The list of them is
sent with every question, and the language model chooses from the descriptions in it.

<figure class="diagram">
<div class="aq-steps" role="group" aria-label="Step to show">
  <button v-for="(b, i) in steps" :key="i" type="button" :class="{ on: step === i }" :aria-pressed="step === i" @click="step = i">{{ b.name }}</button>
</div>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 392" role="img" aria-label="The path of one question: it is sent to the language model with the list of functions, the language model names one, that function is run on the user's machine and reads the saved files, and its result goes back for the answer">
  <rect class="aq-machine" x="14" y="26" width="536" height="352" rx="7" />
  <text class="aq-cap" x="16" y="16">THE USER'S MACHINE</text>
  <text class="aq-cap" x="670" y="16" text-anchor="middle">OVER THE NETWORK</text>
  <g :class="{ 'aq-dim': !lit(1, 2, 3, 4, 5) }">
    <rect class="aq-box" x="34" y="46" width="196" height="150" rx="4" />
    <text class="aq-title" x="132.0" y="117.0" text-anchor="middle">Chat or AI assistant</text>
    <text class="aq-sub" x="132.0" y="136.0" text-anchor="middle">where the question is typed</text>
  </g>
  <g :class="{ 'aq-dim': !lit(1, 2, 4, 5) }">
    <rect class="aq-box aq-llm-box" x="590" y="46" width="160" height="150" rx="4" />
    <text class="aq-title" x="670.0" y="126.0" text-anchor="middle">Language model</text>
  </g>
  <g :class="{ 'aq-dim': !lit(1) }">
    <line class="aq-net" x1="230" y1="72" x2="582" y2="72" />
    <polygon class="aq-net-head" points="590,72 581,67.5 581,76.5" />
    <text class="aq-label" x="410.0" y="65" text-anchor="middle"><tspan class="aq-num">1</tspan>  question, with the list of functions</text>
  </g>
  <g :class="{ 'aq-dim': !lit(2) }">
    <line class="aq-net" x1="590" y1="107" x2="238" y2="107" />
    <polygon class="aq-net-head" points="230,107 239,102.5 239,111.5" />
    <text class="aq-label" x="410.0" y="100" text-anchor="middle"><tspan class="aq-num">2</tspan>  function to call</text>
  </g>
  <g :class="{ 'aq-dim': !lit(4) }">
    <line class="aq-net" x1="230" y1="142" x2="582" y2="142" />
    <polygon class="aq-net-head" points="590,142 581,137.5 581,146.5" />
    <text class="aq-label" x="410.0" y="135" text-anchor="middle"><tspan class="aq-num">4</tspan>  result of the function</text>
  </g>
  <g :class="{ 'aq-dim': !lit(5) }">
    <line class="aq-net" x1="590" y1="177" x2="238" y2="177" />
    <polygon class="aq-net-head" points="230,177 239,172.5 239,181.5" />
    <text class="aq-label" x="410.0" y="170" text-anchor="middle"><tspan class="aq-num">5</tspan>  answer</text>
  </g>
  <g :class="{ 'aq-dim': !lit(1, 2, 3) }">
    <rect class="aq-box" x="34" y="244" width="330" height="120" rx="4" />
    <text class="aq-title" x="48" y="266">Functions</text>
    <text class="aq-sub" x="124" y="266">a name, a description and arguments each</text>
    <rect class="aq-pick" :class="{ on: step === 2 || step === 3 }" x="42" y="276" width="314" height="22" rx="3" />
    <text class="aq-mono" x="50" y="292">explain_cell</text>
    <text class="aq-sub" x="146" y="292">why a cell was assigned to its class</text>
    <text class="aq-mono" x="50" y="314">explain_spot</text>
    <text class="aq-sub" x="146" y="314">why a spot was assigned to its cell</text>
    <text class="aq-mono" x="50" y="336">find_cells</text>
    <text class="aq-sub" x="146" y="336">the cells matching a filter</text>
    <text class="aq-sub" x="50" y="356">and the rest of them</text>
  </g>
  <g :class="{ 'aq-dim': !lit(3) }">
    <line class="aq-local" x1="132" y1="196" x2="132" y2="236" />
    <polygon class="aq-local-head" points="132,244 127.5,235 136.5,235" />
    <text class="aq-label" x="142" y="225"><tspan class="aq-num">3</tspan>  runs the function</text>
    <line class="aq-local" x1="364" y1="304" x2="402" y2="304" />
    <polygon class="aq-local-head" points="410,304 401,299.5 401,308.5" />
    <text class="aq-label" x="386" y="296" text-anchor="middle">reads</text>
    <rect class="aq-box aq-files" x="410" y="281" width="124" height="46" rx="4" />
    <text class="aq-title" x="472.0" y="309.0" text-anchor="middle">Saved files</text>
  </g>
</svg>
<figcaption>{{ steps[step].says }}</figcaption>
</figure>

1. The question is sent to the language model, together with the list of functions it
   may call. An entry of the list is a name, a description and the arguments the
   function takes. This is the entry for `explain_cell`, the description shortened:

   ```
   name         explain_cell
   description  Why a cell was given its class, gene by gene. [...] Use this for
                questions like 'why is cell 2413 Ndnf Gaba' or 'why is cell 18223
                not CA1'.
   arguments    label     The cell label, as in the segmentation.
                vs_class  Class to compare against. Defaults to the runner up.
                top_n     How many genes to list for each side, default 10.
   ```

   An AI assistant receives the same list when it connects to `pciseq-mcp` or to
   pciSeq Viewer.
2. The language model replies with the function to call and its arguments, for example
   `explain_cell` for cell 5016.
3. The function is run on the user's machine. It reads the saved files, or, in the
   live viewer's chat, the fit in progress.
4. The result of the function is sent to the language model.
5. The language model writes the answer from that result.

Steps 2 to 4 are repeated when an answer needs more than one function.

In every way the agent can also show the source code behind a quantity, for example
the adjustment for anisotropy, and in every way the code shown is the code that produced
the numbers. The live viewer's chat reads the installed files, which are the code
running the fit. For a completed run, `pciseq-mcp` and pciSeq Viewer read the file at
the commit the run recorded, from the pciSeq repository on GitHub. This needs a network
connection; the request names a file and a commit and carries nothing from the run. A run that does not record its commit is refused rather than
shown another version. The [code map](./api/code-map.md) names the function and line
for each quantity.

Where an API key is needed, it is saved on the user's machine. The files `pciSeq.fit`
saved stay on that machine too; only the question and the results of the
functions are sent to the language model.
