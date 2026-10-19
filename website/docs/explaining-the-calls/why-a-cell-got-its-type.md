
# How a cell's call was made

A cell is assigned the class with the highest posterior probability. For every class the
model adds three terms, and a softmax over all classes turns the totals into
probabilities:

$$
\text{score}_k = \underbrace{\text{gene log-likelihood}_k}_{\text{the cell's counts}}
+ \underbrace{\log \pi_k}_{\text{the prior}}
+ \underbrace{\text{MRF}_k}_{\text{the neighbours}}
$$

[`check_cell`](../api/reference.md#check-cell) compares the assigned class of a cell with
a second class of your choice and shows how each of these terms contributes. It reads the
values the model used in its last cell class update, so the figure agrees with the
probabilities in `cellData`.

The examples on this page use the 3D dataset described in the [overview](overview.md),
fitted with `mrf_beta = 1.5`. The second half of the page compares the same cell in a
fit without the spatial term, `mrf_beta = 0`.

## The cell

Cell 18223 lies in the dentate gyrus, at the tip of the upper blade. Panel **a** is the
DAPI image of the section, with the cell circled. Panel **b** is a close-up of the same
cell, with its outline in red and, in blue, the nine nearest cells that the spatial term
reads. Neighbours are taken in 3D, so one of the nine sits on another plane and is not
drawn.

<figure class="diagram" id="fig-the-cell">
<img src="/explaining-the-calls/cell-18223-map.png" alt="Cell 18223 in the section and up close">
<figcaption><strong>Figure 3.1.</strong> Cell 18223 in the section, panel <strong>a</strong>,
and up close with its nine nearest cells, panel <strong>b</strong>.</figcaption>
</figure>

Without the spatial term the cell is assigned `030 L6 CT CTX Glut`, a cortical class. With
it, the cell is assigned `037 DG Glut`. The rest of this page examines the score of both
fits, term by term and iteration by iteration.

## Calling `check_cell`

Cell 18223 has a total gene count of 39 and is assigned `037 DG Glut` with probability
1.00. The following compares it with `030 L6 CT CTX Glut`, a cortical class:

```python
import pandas as pd

# the fitted model, saved by pciSeq.fit at <output_path>/pciSeq/data/debug/pciSeq.pickle
obj = pd.read_pickle('pciSeq.pickle')
obj.check_cell(18223, '030 L6 CT CTX Glut')
```

`fit` saves the fitted model when `save_data` is `True`, the default, see
[the fitted model](overview.md#the-fitted-model).
The cell is given by its label in the segmentation passed to `fit`.

<figure class="diagram" id="fig-check-cell-mrf">
<img src="/explaining-the-calls/cell-18223-mrf.png" alt="check_cell output for cell 18223">
<figcaption><strong>Figure 3.2.</strong> <code>check_cell</code> for cell 18223 in the fit
with the spatial term.</figcaption>
</figure>

The figure has four panels.

1. **Top left.** For each gene, the log-likelihood of the cell's count under the assigned
   class minus its log-likelihood under the second class. A positive difference means the
   count is more likely under the assigned class, so the gene favours it. The panel shows
   up to `top_n` genes with the largest positive differences, and the title gives their
   sum.
2. **Top right.** The same difference for the genes where it is negative: the count is
   more likely under the second class, so the gene favours the second class. The panel
   shows up to `top_n` genes with the most negative differences.
3. **Bottom left.** The three terms of the score for both classes: the gene log-likelihood
   summed over all genes (the top panels break this term down by gene), the log prior and
   the MRF term. Higher is better: the class with
   the larger total wins.
4. **Bottom right.** The posterior over all classes, for the `top_classes` most likely
   classes and the second class. The note gives the probability of the classes not shown.

For cell 18223:

- **Top panels.** Synpr, Sema5a and Pde1a favour `037 DG Glut`, by 2.38, 1.47 and 0.94;
  the ten genes on that side sum to 8.26. Neurod6 and Rgs4 favour `030 L6 CT CTX Glut`, by
  1.75 and 1.38; the ten genes on that side sum to -4.99.
- **Bottom left.** The gene log-likelihood is -99.7 for `037 DG Glut` against -104.4 for
  `030 L6 CT CTX Glut`, so the genes alone favour DG by 4.6. The priors are equal. The MRF
  term is 13.5 for `037 DG Glut` and 0 for `030 L6 CT CTX Glut`: all nine neighbours of the
  cell are DG cells, and 13.5 is the largest value the term takes with nine neighbours.
- **Bottom right.** The posterior is 100% on `037 DG Glut`.

Both terms favour `037 DG Glut`, and the posterior is 1.00.

## The returned table

```python
table, contr, fig = obj.check_cell(18223, '030 L6 CT CTX Glut')
```

`check_cell` returns three things: a table with one row for each gene in the two gene charts
of the figure (top left and top right), the per-gene log-likelihoods the charts are drawn
from, and the figure itself. The first is shown here:

<!--@include: ./_tables/cell-18223-mrf.md-->

Its columns are:

- **Cells typed as a class, mean counts.** Estimated from the spots, across the run: the average counts of
  the gene over all cells, each weighted by its probability of the class.
- **Predicted counts for cell 18223.** Expected, for this cell: the counts of the gene the
  model expects if the cell were of that class, the cell type definition rescaled by the
  [scaling factors](../how-it-works/warping-the-reference.md), see the
  [example below](#the-predicted-counts-are-not-the-mean-counts). A gene favours the class whose
  prediction is closer to those counts.
- **Cell 18223, counts.** Estimated from the spots, in this cell: the counts of the gene.

Synpr is the strongest gene. The Synpr count of the cell is 1.7. `037 DG Glut` predicts
0.76 and `030 L6 CT CTX Glut` predicts 0.10. The count is small, but the prediction under
`030 L6 CT CTX Glut` is near zero, so each Synpr spot adds 1.74 in favour of
`037 DG Glut`. With a count of 1.7 Synpr contributes 2.38, the largest bar in the
top-left panel. The cells assigned to `037 DG Glut` have a mean Synpr count of 0.61,
those assigned to `030 L6 CT CTX Glut` 0.05.

::: details How the 2.38 is computed
The gene log-likelihood is a negative binomial with mean $\mu$, the expected count, and
dispersion $r$, the `rSpot` setting, here 2. For a count $x$:

$$
\log \text{NB}(x;\, r, \mu) = x \log\frac{\mu}{r+\mu} + r \log\frac{r}{r+\mu} + \log\frac{\Gamma(x+r)}{x!\,\Gamma(r)}
$$

The last term does not depend on $\mu$ and cancels in the difference between two classes
$A$ and $B$:

$$
\Delta = x \left[\log\frac{\mu_A}{r+\mu_A} - \log\frac{\mu_B}{r+\mu_B}\right] + r \log\frac{r+\mu_B}{r+\mu_A}
$$

The first term grows with the count, the second does not. For Synpr, with the expected
counts from the table:

| class | $\mu$ | $\mu/(r+\mu)$ | $\log$ |
| --- | --- | --- | --- |
| `037 DG Glut` | 0.765 | 0.2767 | -1.2850 |
| `030 L6 CT CTX Glut` | 0.102 | 0.0487 | -3.0212 |

Each spot adds $-1.2850 - (-3.0212) = 1.7362$ in favour of `037 DG Glut`. The second term
is $2 \log(2.102/2.765) = -0.548$ and favours `030 L6 CT CTX Glut`, which predicts fewer
spots overall. With $x = 1.689$:

$$
\Delta = 1.689 \times 1.7362 - 0.548 = 2.38
$$

With a Synpr count of 0, $\Delta = -0.55$: the absence of a gene favours the class that
predicts fewer.
:::

### The predicted counts are not the mean counts

For Synpr the model predicts a count of 0.76 in cell 18223 as `037 DG Glut`, while the cells
typed as `037 DG Glut` have 0.61 on average. The two numbers are different quantities.

The prediction is built from the single-cell reference:

| factor | value |
| --- | --- |
| reference mean of Synpr in `037 DG Glut` | 24.68 |
| × `Inefficiency` | 0.1 |
| × eta of Synpr | 0.52 |
| × theta of cell 18223 as `037 DG Glut` | 0.52 |
| + `SpotReg` | 0.1 |
| **prediction** | **0.76** |

gamma, the factor for a given gene in a given cell, is not part of the product. It expresses the
discrepancy between the count in the cell and the prediction, regularised by
`rSpot`.

The prediction and the mean counts need not agree. eta is one number per gene, shared by
all classes, so it cannot correct the reference for each class separately, and the
reference does not match the in situ data exactly. For Sema5a in the same table the
prediction as `037 DG Glut` is 3.96 against mean counts of 6.94. What decides the
call is the comparison within the cell: the observed Synpr count of 1.69 is closer to the
prediction under `037 DG Glut`, 0.76, than under `030 L6 CT CTX Glut`, 0.10.

## The same cell without the spatial term

Fitted with `mrf_beta = 0`, and otherwise identical settings, cell 18223 is assigned
`030 L6 CT CTX Glut`:

```python
# the fitted model of the run with mrf_beta = 0
obj_nomrf = pd.read_pickle('pciSeq_nomrf.pickle')
obj_nomrf.check_cell(18223, '037 DG Glut')
```

<figure class="diagram" id="fig-check-cell-nomrf">
<img src="/explaining-the-calls/cell-18223-nomrf.png" alt="check_cell output for cell 18223 without the MRF">
<figcaption><strong>Figure 3.3.</strong> The same cell in the fit without the spatial
term.</figcaption>
</figure>

- **Top panels.** Neurod6, Rgs4 and Rprm favour `030 L6 CT CTX Glut` by 5.27, 3.54 and
  1.24. Pde1a, Trp53i11 and Sema5a favour `037 DG Glut` by 0.74, 0.66 and 0.48.
- **Bottom left.** The gene log-likelihood is -95.2 for `030 L6 CT CTX Glut` against
  -107.6 for `037 DG Glut`: the genes favour L6 CT by 12.4. There is no MRF term. Seven of
  the nine neighbours are DG cells in this fit, but they do not enter the score.
- **Bottom right.** The posterior is 98.7% `030 L6 CT CTX Glut`; `037 DG Glut` is below
  0.1%.

In both fits the gene log-likelihood is larger for the assigned class. What differs is not how the evidence
is weighed but the evidence itself: the two fits give the cell different spots.

| gene | favours | gene counts, without MRF | gene counts, with MRF |
| --- | --- | --- | --- |
| Neurod6 | L6 CT | 3.49 | 1.41 |
| Rgs4 | L6 CT | 4.13 | 1.92 |
| Rprm | L6 CT | 1.53 | 0.82 |
| Synpr | DG | 0.11 | 1.69 |
| Sema5a | DG | 3.91 | 6.51 |

The total count is almost the same, 38 against 39. Scoring the spots the cell
holds in the fit with the MRF under the parameters of the fit without it turns the
preference of the genes from 12.4 for L6 CT to 4.4 for DG. The spots the cell holds
account for almost all of the difference.

### How the two fits diverge

`check_cell` reads the last iteration only. The table below follows cell 18223 through
the iterations of both fits. Class is the most likely class of the cell at that iteration
and Prob its posterior probability.

<table>
<thead>
<tr><th rowspan="2">iteration</th><th colspan="3" style="text-align: center">without MRF</th><th colspan="4" style="text-align: center">with MRF</th></tr>
<tr><th>class</th><th>Prob</th><th>Δ log-lik¹</th><th>class</th><th>Prob</th><th>Δ log-lik¹</th><th>neighbours²</th></tr>
</thead>
<tbody>
<tr><td>0</td><td>DG</td><td>0.77</td><td>+1.2</td><td>DG</td><td>0.77</td><td>+1.2</td><td>8 DG, 1 L6 CT</td></tr>
<tr><td>1</td><td>L6 CT</td><td>0.62</td><td>-0.5</td><td>DG</td><td>1.00</td><td>-0.5</td><td>9 DG</td></tr>
<tr><td>2</td><td>L6 CT</td><td>0.98</td><td>-6.2</td><td>DG</td><td>1.00</td><td>+1.8</td><td>9 DG</td></tr>
<tr><td>10</td><td>L6 CT</td><td>0.97</td><td>-10.6</td><td>DG</td><td>1.00</td><td>+4.4</td><td>9 DG</td></tr>
<tr><td>last</td><td>L6 CT</td><td>0.99</td><td>-12.4</td><td>DG</td><td>1.00</td><td>+4.6</td><td>9 DG</td></tr>
</tbody>
<tfoot>
<tr><td colspan="8"><div style="width: 0; min-width: 100%; font-size: 0.92em">¹ Gene log-likelihood of <code>037 DG Glut</code> minus that of <code>030 L6 CT CTX Glut</code>: positive values mean the gene counts alone, without the prior and the MRF term, favour DG; negative values favour L6 CT.<br>² Most likely class of the nine nearest neighbours of the cell at the end of the iteration.</div></td></tr>
</tfoot>
</table>

At the start of iteration 0 each spot is split equally between its nine nearest cells and
the background, regardless of distance or gene. The call in iteration 0 is made from these
counts and is the same in both fits. At the end of iteration 0 the spots are reassigned,
this time using distance and expression. The Neurod6 count of the cell goes from 0.6 to
1.9 and the Rgs4 count from 1.7 to 2.6, two genes that favour L6 CT. At iteration 1 the genes
therefore favour L6 CT by 0.5, in both fits.

Without the MRF term the gene log-likelihood is the only term that differs, and the cell
is assigned L6 CT. With it, the MRF term, from eight DG neighbours out of nine, outweighs
the difference of 0.5 and the cell stays DG. From then on the
counts of each fit follow its class: as L6 CT the Neurod6 and Rgs4 counts rise and the
Synpr and Sema5a counts fall, as DG the reverse. In the fit with the
MRF the MRF term changes the class at iteration 1 only; from iteration 2 the gene
log-likelihood favours DG on its own.
