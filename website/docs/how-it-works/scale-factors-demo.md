---
# standard column, margins and reading measure. Only the demo runs wide.
pageClass: demo-wide
---

# 2a. Demo: the scaling factors

<DemoFrame src="/cells-demo.html?narrow=1&amp;fold=1" wide
           title="Interactive demo of the scaling factors" />

**Simplifications.**

**Eta** is held at 1 throughout, so the demo shows three factors rather than four. In a run,
`rGene`, 20 by default, sets how strongly the prior holds it at 1.

The class prior and the spatial term are both left out, which means the
probabilities come from the likelihood alone; the
[cell-class assignment](../the-model/cell-class.md) page gives the full expression.
