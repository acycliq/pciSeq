# 4. Assigning spots to cells

Each spot is assigned to one of its nearest cells or to the background, as a
probability. This is the last step of the sweep: the new assignments change the gene
counts per cell, which are the input to the next sweep.

## The assignment score

For a spot and a candidate cell, at the default settings, the score is a sum of five
terms. One depends on the spot's position, four on the expected count of its gene in the
cell. A sixth term, the inside-cell bonus, adds the value of `InsideCellBonus` when the
spot falls inside the candidate cell's segmentation mask. It is 0 by default and is
omitted from the diagram and from the rest of this page, see
[configuration](../api/configuration.md#insidecellbonus).

<figure class="diagram">
<svg class="sb-svg" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 940 338" role="img" aria-label="The terms of the spot-to-cell score">
<rect class="sb-root" x="380" y="12" width="180" height="42" rx="11" />
<text class="sb-root-txt" x="470" y="39" text-anchor="middle">spot <tspan font-size="19">&#8594;</tspan> cell score</text>
<path class="sb-link" d="M 470,54 L 470,92" />
<path class="sb-link" d="M 122,92 L 586,92" />
<path class="sb-link" d="M 122,92 L 122,110" />
<path class="sb-link" d="M 586,92 L 586,110" />
<rect class="sb-panel-where" x="24" y="110" width="196" height="216" rx="12" />
<text class="sb-head-where" x="122" y="134" text-anchor="middle">POSITION</text>
<text class="sb-sub" x="122" y="150" text-anchor="middle">1 term</text>
<rect class="sb-panel-what" x="256" y="110" width="660" height="216" rx="12" />
<text class="sb-head-what" x="586" y="134" text-anchor="middle">EXPECTED COUNT</text>
<text class="sb-sub" x="586" y="150" text-anchor="middle">4 terms</text>
<text class="sb-plus" x="238" y="247" text-anchor="middle">+</text>
<text class="sb-plus" x="426" y="247" text-anchor="middle">+</text>
<text class="sb-plus" x="586" y="247" text-anchor="middle">+</text>
<text class="sb-plus" x="746" y="247" text-anchor="middle">+</text>
<rect class="sb-card sb-card-where" x="44" y="174" width="156" height="136" rx="10" />
<circle class="sb-badge-where" cx="122" cy="208" r="21" />
<text class="sb-badge-txt" x="122" y="213" text-anchor="middle" font-size="12">loglik</text>
<text class="sb-name" x="122" y="252" text-anchor="middle">Gaussian fit</text>
<text class="sb-q" x="122" y="275" text-anchor="middle">distance to the</text>
<text class="sb-q" x="122" y="292" text-anchor="middle">centroid</text>
<rect class="sb-card sb-card-what" x="272" y="174" width="148" height="136" rx="10" />
<rect class="sb-badge-what" x="311" y="190" width="70" height="36" rx="18" />
<text class="sb-badge-txt" x="346" y="214" text-anchor="middle" font-size="18"><tspan font-size="12">from </tspan>&#956;</text>
<text class="sb-name" x="346" y="252" text-anchor="middle">Class expression</text>
<text class="sb-q" x="346" y="275" text-anchor="middle">expected count</text>
<text class="sb-q" x="346" y="292" text-anchor="middle">of the gene</text>
<rect class="sb-card sb-card-what" x="432" y="174" width="148" height="136" rx="10" />
<rect class="sb-badge-what" x="471" y="190" width="70" height="36" rx="18" />
<text class="sb-badge-txt" x="506" y="214" text-anchor="middle" font-size="18"><tspan font-size="12">from </tspan>&#952;</text>
<text class="sb-name" x="506" y="252" text-anchor="middle">Cell scale</text>
<text class="sb-q" x="506" y="275" text-anchor="middle">counts of the cell</text>
<text class="sb-q" x="506" y="292" text-anchor="middle">vs its class</text>
<rect class="sb-card sb-card-what" x="592" y="174" width="148" height="136" rx="10" />
<rect class="sb-badge-what" x="631" y="190" width="70" height="36" rx="18" />
<text class="sb-badge-txt" x="666" y="214" text-anchor="middle" font-size="18"><tspan font-size="12">from </tspan>&#947;</text>
<text class="sb-name" x="666" y="252" text-anchor="middle">Cell-gene scale</text>
<text class="sb-q" x="666" y="275" text-anchor="middle">this gene in the</text>
<text class="sb-q" x="666" y="292" text-anchor="middle">cell vs its class</text>
<rect class="sb-card sb-card-what" x="752" y="174" width="148" height="136" rx="10" />
<rect class="sb-badge-what" x="791" y="190" width="70" height="36" rx="18" />
<text class="sb-badge-txt" x="826" y="214" text-anchor="middle" font-size="18"><tspan font-size="12">from </tspan>&#951;</text>
<text class="sb-name" x="826" y="252" text-anchor="middle">Gene efficiency</text>
<text class="sb-q" x="826" y="275" text-anchor="middle">detection rate</text>
<text class="sb-q" x="826" y="292" text-anchor="middle">of the gene</text>
</svg>
<figcaption>The terms of the score. One depends on the spot's position, four on the expected count of its gene in the cell. They are added.</figcaption>
</figure>

### Gaussian fit

Each cell is described by a Gaussian centred on its fitted centroid, with the cell's
fitted covariance. The term is the log density of the spot's position under that
Gaussian: it is highest at the centre and decreases with distance, more rapidly along
the cell's short axes. It does not depend on the gene.

### Expression fit

Together, the four expression terms are the log of the expected count of the spot's gene
in the candidate cell, taken from the [warped definitions](warping-the-reference.md);
each term is derived from one of its factors. The first three depend on the cell's class,
which is not known exactly, so each is averaged over the cell's class probabilities from
[cell typing](cell-to-celltype.md). The fourth depends on the gene alone.

- **Class expression.** Derived from mu, the expected count of the gene in each class:
  its log, averaged over the cell's class probabilities. It is high only when the
  probabilities are concentrated on classes that express the gene. Formally, it is the
  dot product of two vectors over the classes, the cell's class probabilities and the
  log expected counts of the gene, and it measures how closely the two coincide. This
  term is therefore referred to here as the alignment. Two cells with the same class
  probabilities have the same alignment. When the cell's likely class expresses the gene
  more than the alternatives, a cell less confident of that class has a lower alignment,
  because the probability it places on classes with little of the gene enters through
  the log.
- **Cell scale.** Derived from theta, the cell's total observed count relative to what
  its class predicts: its log, averaged in the same way. It does not depend on the gene:
  a cell holding more transcripts than its class predicts scores higher for every gene,
  and so attracts nearby spots of any gene. This term is referred to here as the
  gravity. Unlike physical gravity, it does not depend on distance; distance enters only
  through the Gaussian fit.
- **Cell-gene scale.** Derived from gamma, this gene's observed count in the cell
  relative to what its class and theta predict: its log, averaged in the same way. When
  gamma is above 1, the cell holds more of the gene than predicted; below 1, less. This
  term is referred to here as the enrichment, a name that covers depletion as well.
  Between two cells with the same class probabilities and the same theta, it is the only
  expression term that differs.
- **Gene efficiency.** Derived from eta, the gene's detection efficiency: its log. It is
  the same for every candidate cell, so it does not choose between cells. It enters the
  comparison with the background, which has no efficiency term, and lowers the score of
  a poorly detected gene against it.

### An intuitive reading of the terms

The role of each term can be seen without the formula. Suppose a spot lies between two
candidate cells and every term but one is equal for the two. The remaining term decides,
and in each case it gives the answer one would give by inspection. For simplicity, the
shape and the centroid of each cell are taken as fixed.

- **Gaussian fit.** Take two cells of the same shape that differ only in position. The
  spot goes to the nearer cell.
- **Alignment.** Consider a Plp1 spot equidistant from two cells of the same shape and
  total count, one an oligodendrocyte and the other a CA1 pyramidal cell. Plp1 is an
  oligodendrocyte gene, so the spot goes to the oligodendrocyte. If both cells are
  oligodendrocytes, it goes to the cell more confidently typed as one.
- **Gravity.** Place a spot at the same distance from two oligodendrocytes of the same
  shape, equally confidently typed, one holding more transcripts in total than the
  other. The spot goes to the cell with more transcripts: a cell that has produced more
  reads is the more likely source of one more. Because the two cells share a class, more
  transcripts means a higher theta. Between cells of different classes, theta compares
  each cell with what its own class predicts, not with the other cell's count.
- **Enrichment.** Suppose a Plp1 spot is equally close to two oligodendrocytes alike in
  shape, typing confidence and total count, one holding more Plp1 than the other. The
  spot goes to the cell already holding more Plp1.
- **Gene efficiency.** Eta is the same for both cells and cannot decide between them. It
  acts on the comparison with the background: a gene that is poorly detected yields few
  reads from cells, so a larger share of its reads are misreads, and a read of it is
  more readily assigned to the background.

For a real spot, several terms differ at once, and their sum decides. The tie isolates
one term to show its direction.

## The background option

The background is scored with the gene's [misread density](misread-density.md), the
same at every position. The more a candidate cell's score falls below the background's,
the smaller its share of the spot's probability; for a spot that no candidate cell
explains better than the background, the background receives the largest share. This is
how misreads are removed.

## Soft assignments

A softmax over the candidates and the background gives the spot a probability for each,
for example 0.9 on one cell and 0.1 on a neighbour. The gene counts for the next sweep
are accumulated with these probabilities, so an ambiguous spot is shared between cells
rather than committed to one.

## Feedback into the next sweep

The updated counts per cell are the input to the [misread density](misread-density.md)
and the [warping](warping-the-reference.md) of the next sweep. The loop stops when the
spot probabilities stop changing, see [convergence](overview.md#convergence).
