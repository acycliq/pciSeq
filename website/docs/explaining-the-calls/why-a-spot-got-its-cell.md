
# How a spot's call was made

A spot is assigned to one of its `nNeighbors` nearest cells or to the background, as a
probability. Each candidate cell gets a score, the background gets one too, and a softmax
over them gives the probabilities:

$$
\begin{aligned}
\text{score}_c = \;& \text{Gaussian fit} && \text{position} \\
+\;& \text{class expression} + \text{cell scale} + \text{cell-gene scale} + \text{gene efficiency} && \text{expression} \\
+\;& \text{inside-cell bonus} && \text{segmentation}
\end{aligned}
$$

The background scores the [misread density](../how-it-works/misread-density.md) of the
spot's gene, the same at every position. The terms are described in
[assigning spots to cells](../how-it-works/spots-to-cells.md).

[`check_spot`](../api/reference.md#check-spot) shows these terms for every candidate of a
spot. It reads the values the model used in its last spot update, so what it shows agrees
with the probabilities in `geneData`.

The examples use the same two fits as the [cell page](why-a-cell-got-its-type.md), with
`mrf_beta = 1.5` and `mrf_beta = 0`, and one spot of cell 18223.

## Calling `check_spot`

The example uses the same cell and the same two fits as
[Figure 3.1](why-a-cell-got-its-type.md#fig-the-cell). The question is now on the spot
side: which cell does the spot belong to, and how the answer changes when the spatial
term (mrf) is switched off.

Spot 1642419 is a Synpr spot. In the fit with the spatial term it goes to cell 18223:

```python
import pandas as pd

obj = pd.read_pickle('pciSeq.pickle')
obj.check_spot(1642419)
```

The spot is given by its index in the spots table passed to `fit`.

<figure class="diagram" id="fig-check-spot-scores">
<img src="/explaining-the-calls/spot-1642419-mrf-scores.png" alt="check_spot score decomposition for spot 1642419">
<figcaption><strong>Figure 3.4.</strong> The score of every candidate of spot 1642419,
term by term.</figcaption>
</figure>

The bars are the score of each candidate, term by term, with the candidates ordered by
distance and the background last. The black tick is the total. A taller total means a
better explanation of the spot.

<figure class="diagram" id="fig-check-spot-probs">
<img src="/explaining-the-calls/spot-1642419-mrf-probs.png" alt="check_spot assignment probabilities for spot 1642419">
<figcaption><strong>Figure 3.5.</strong> The assignment probabilities of the same spot,
the softmax of the totals in Figure 3.4.</figcaption>
</figure>

The second chart is the softmax of those totals, which is the assignment the model uses:
0.74 on cell 18223, 0.13 on cell 21574, 0.10 on cell 17371 and 0.01 on the background.

Each expression term is a log, so its sign says whether the quantity inside the log is
above or below 1.

Take the class expression, the orange segment of every bar in the score chart above. It sits
above zero, so the quantity inside its log is above 1. The term is
`log(mean expression x Inefficiency + SpotReg)`, averaged over the cell's class
probabilities. Every candidate here is `037 DG Glut` with probability near 1, so the
average reduces to the value for DG. Synpr has a mean of 24.68 in `037 DG Glut`, and with
`Inefficiency` 0.1 and `SpotReg` 0.1 that gives `log(24.68 x 0.1 + 0.1) = log(2.57) = 0.943`, the value every DG
candidate shows here.

It is not positive by construction: the sign follows the class mean. For a cell
confidently of one class, as here, the term is zero when
`mean expression x Inefficiency + SpotReg` equals 1, which in this run is a class mean of
9 counts, and negative below that. For a cell spread over several classes it is the
average of those logs.

The cell scale and the cell-gene scale are positive when the cell holds more than its class
predicts. Here the cell-gene scale is positive for cell 18223, 0.18, and cell 21574, 0.07, the
cells with more Synpr than their class predicts; the cell scale is positive for cells 22339,
0.28, and 21574, 0.10, which hold more reads in total than their class predicts.

## The returned table

`check_spot` returns the same numbers as a table, one row per candidate and a final row
for the background:

<!--@include: ./_tables/spot-1642419-mrf.md-->

<p class="table-note" id="table-3-1"><strong>Table 3.1.</strong> The output of
<code>check_spot</code> for spot 1642419, in the fit with the spatial term.</p>

- **Gaussian fit**. The log density of the spot's position under the cell's Gaussian
  distribution. For simplicity, the fits of this example were run without the
  [cell geometry](../the-model/gaussian.md) update, so the Gaussians are fixed: spherical,
  with the same covariance for every cell and centred on the segmentation centroid. Only
  the distance to the cell's centroid counts then, and the three axes are independent. With
  the update, the centre and the covariance follow the spots and the term changes with
  them.
- **class expression, the alignment**. How well the cell's likely class matches the gene. When several cells
  compete for the same spot and everything else is equal, it goes to the cell whose likely
  class expresses that gene most. Two cells with the same class probabilities get the same
  value, since the gene is the same for every candidate of a spot, which is why every
  candidate here shows 0.943.
- **cell scale, the gravity**. The factor behind it is theta: how far a cell's total
  counts sit above or below what its likely class predicts. When several cells compete for
  the same spot and everything else is equal, it goes to the cell already holding more
  reads in total. Only the totals count here, so the term is the same whatever the gene of
  the spot.
- **cell-gene scale, the enrichment**. The factor behind it is gamma: how much of a gene
  is observed in a cell relative to the amount expected for its likely class. When several
  cells compete for the same spot and all other factors are equal, the spot is assigned to
  the cell with the higher observed-to-expected ratio for that gene.
- **gene efficiency**. The factor behind it is eta: how well the gene is detected across
  the whole experiment. It is the same for every candidate cell, so it cannot choose
  between them. It does weigh the cells against the background, which carries no such
  term: a poorly detected gene lowers every cell's score and makes the background more
  likely.
- **misread**. The background's whole score: how likely a spot of this gene is to be a
  misread. It is the same at every position, so a cell has to beat it or the spot goes to
  the background.
- **sum, prob**. The total and its softmax over the rows.

<p class="table-note">Gene efficiency and misread density are not the same thing. The
first scales what a cell is expected to hold of the gene, so it sits in every cell's
score. The second is the background's own rate for that gene, learned from the spots no
cell explains. A gene can be well detected and still produce many misreads.</p>

<p class="table-note">The inside-cell bonus is not shown: <code>InsideCellBonus</code> is 0
in these runs, so the column was empty. When it is set, a spot whose pixel falls inside the
cell's segmentation gets that bonus added to its score.</p>

Cell 18223 is the nearest candidate: its Gaussian fit is 1.7 above cell 17371 and 2.4
above cell 21574. It does not win on the expression terms. It has the largest cell-gene
scale, 0.18 against 0.07 for cell 21574 and negative values for the rest, but a lower
cell scale, -0.66 against 0.10, and over the three expression terms together cell 21574
is ahead, 1.11 against 0.47. Distance carries the call.

## The same spot without the spatial term

In the fit with `mrf_beta = 0`, cell 18223 is `030 L6 CT CTX Glut` rather than
`037 DG Glut`. The spot has not moved and the Gaussians are fixed, so its position terms
are unchanged:

```python
obj_nomrf = pd.read_pickle('pciSeq_nomrf.pickle')
obj_nomrf.check_spot(1642419)
```

<!--@include: ./_tables/spot-1642419-nomrf.md-->

<p class="table-note" id="table-3-2"><strong>Table 3.2.</strong> The same spot in the
fit without the spatial term. The Gaussian fit column is identical to
<a href="#table-3-1">Table 3.1</a>.</p>

<div class="two-tables">
<div>
<table>
<thead><tr><th>term for cell 18223</th><th>without MRF</th><th>with MRF</th></tr></thead>
<tbody>
<tr><td>Gaussian fit</td><td>-10.405</td><td>-10.405</td></tr>
<tr><td>class expression</td><td>-2.203</td><td>+0.943</td></tr>
<tr><td>sum</td><td>-14.149</td><td>-10.590</td></tr>
<tr><td>prob</td><td>0.06</td><td>0.74</td></tr>
</tbody>
</table>
<p class="table-note" id="table-3-3"><strong>Table 3.3.</strong> The terms of cell 18223
that change between the two fits.</p>
</div>
<div>
<table>
<thead><tr><th>cell</th><th>Cplx2</th><th>Sema5a</th><th>Snca</th><th>Nrn1</th><th>Bcl11b</th><th>...</th><th>Synpr</th><th>total</th></tr></thead>
<tbody>
<tr><td>17371</td><td>3.3</td><td>10.8</td><td>3.8</td><td>0.5</td><td>4.3</td><td></td><td>0.5</td><td>60.0</td></tr>
<tr><td>21574</td><td>13.0</td><td>10.9</td><td>4.7</td><td>4.4</td><td>2.3</td><td></td><td>3.0</td><td>98.8</td></tr>
</tbody>
</table>
<p class="table-note" id="table-3-4"><strong>Table 3.4.</strong> What the two other
<code>037 DG Glut</code> candidates hold in the fit without the spatial term. The five
genes with the largest counts, then Synpr, then the cell total.</p>
</div>
</div>

- **`Gaussian fit`**, -10.405 in both fits. With the Gaussians fixed, it depends only on
  where the spot is.
- **`class expression`**, +0.943 with the mrf, -2.203 without it. The term is high when the
  cell's likely class expresses the gene strongly, and the cell is confidently of that
  class. Cell 18223 is now classified as `030 L6 CT CTX Glut`, which does not express
  Synpr, while `037 DG Glut` does.

The spot is now assigned to cell 18223 with probability 6%, against 73.7% with the mrf,
and its most likely parent is cell 21574, a `037 DG Glut` cell, with 49.2%. Cell 21574 is
only the third closest to the spot. Cell 17371 is the second closest and is `037 DG Glut`
too, yet the spot does not go to it, because it holds fewer reads in total, 60.0 against
98.8, and expresses fewer Synpr, 0.5 against 3.0, see [Table 3.4](#table-3-4). That is
intuitively sensible: a spot is more likely to come from a cell that is bigger in
terms of total reads and expresses more of the gene than from one that does not.


Cell 18223 here is the same cell as in
[Figure 3.1](why-a-cell-got-its-type.md#fig-the-cell), the
same fits and the same plane, 57. There it was panel **b**, outlined in red among its
neighbours.

<figure class="diagram" id="fig-spot-panels">
<div class="two-panel">
  <div>
    <img src="/explaining-the-calls/spot-1642419-mrf-map.png" alt="the three candidate cells with the spatial term">
    <p>with the spatial term</p>
  </div>
  <div>
    <img src="/explaining-the-calls/spot-1642419-nomrf-map.png" alt="the same cells without the spatial term">
    <p>without the spatial term</p>
  </div>
</div>
<figcaption><strong>Figure 3.6.</strong> The same field at plane 57, with Synpr the only
gene shown. Cells carry the colour of their class, and the line joins the spot to the
cell it was assigned to. Cell 18223 is <code>037 DG Glut</code> on the left and <code>030 L6 CT CTX Glut</code>, in
green, on the right, while cells 17371 and 21574 stay <code>037 DG Glut</code> in both.
</figcaption>
</figure>
