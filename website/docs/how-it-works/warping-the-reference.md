
# 2. Warping the cell type definitions

The cell type definitions give the mean expression of each gene in each class. The
observed counts are on a different scale: detection efficiency differs between genes,
transcript yield differs between cells, and the cell type definitions and the experiment
differ in overall scale.
pciSeq rescales the expected expression to the scale of the current experiment before
comparing. The rescaling is the warp.

The step reads the current gene counts per cell, the current cell class estimates and the
raw cell type definitions, and produces a warped expected expression rescaled at every
level. [Cell typing](cell-to-celltype.md) scores cells against classes using it.

The correction factors are nuisance parameters. They are estimated because the cell
classes and spot assignments depend on them, and nothing
observed corresponds to them directly. They are identified through the agreement they
produce between cells and their assigned classes.

## The scaling factors

The warp is a family of four factors, each rescaling the expected expression at a
different level of detail. pciSeq calls them inefficiencies. From the broadest to the
most specific:

- **Inefficiency.** One constant for the whole reference, the `Inefficiency` setting.
  It rescales every gene in every cell type definition to the overall scale of the
  observed counts. It can be below or above 1.

- **eta** ($\eta_g$). One factor per gene, shared by all cells: the gene's detection
  efficiency relative to the constant above.

- **theta** ($\theta_{c\mid k}$). One factor per cell and candidate class: the cell's
  total yield relative to what the class predicts, applied to all its genes.

- **gamma** ($\gamma_{g,c\mid k}$). One factor per gene, cell and candidate class: the
  residual mismatch of a given gene in a given cell that the broader factors leave.

Inefficiency and eta do not depend on the cell's class. theta and gamma do, since the
expectation they correct is class-specific, and they are computed for every candidate
class. theta is what lets [cell typing](cell-to-celltype.md) score a cell against every
class; gamma does not enter cell typing directly: it is integrated out, which gives the
negative binomial.

## Granularity of the factors

Ordered by how much of the experiment each covers: Inefficiency applies to everything,
gamma to a given gene in a given cell under a given class.

<figure class="diagram">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 920 460" role="img" aria-label="A pyramid of the four scaling factors">
  <defs></defs>
  <line x1="70" y1="396" x2="70" y2="66" stroke="currentColor" stroke-opacity="0.5" stroke-width="1.5" />
  <path d="M 70,54 L 76,68 L 70,64 L 64,68 Z" fill="currentColor" fill-opacity="0.5" />
  <text class="ip-axis" transform="rotate(-90 90,230)" x="90" y="230" text-anchor="middle">GRANULARITY INCREASES</text>
  <path class="ip-tier" fill="#34d399" d="M 300,60.0 L 342.9,141.0 L 257.1,141.0 Z" />
  <path class="ip-tier" fill="#10b981" d="M 252.9,149.0 L 347.1,149.0 L 387.9,226.0 L 212.1,226.0 Z" />
  <path class="ip-tier" fill="#059669" d="M 207.9,234.0 L 392.1,234.0 L 432.9,311.0 L 167.1,311.0 Z" />
  <path class="ip-tier" fill="#047857" d="M 162.9,319.0 L 437.1,319.0 L 477.9,396.0 L 122.1,396.0 Z" />
  <line class="ip-leader" x1="343" y1="100.5" x2="540" y2="100.5" />
  <line class="ip-leader" x1="368" y1="187.5" x2="540" y2="187.5" />
  <line class="ip-leader" x1="413" y1="272.5" x2="540" y2="272.5" />
  <line class="ip-leader" x1="458" y1="357.5" x2="540" y2="357.5" />
  <g>
    <text x="550" y="96" class="ip-glyph">&#947; <tspan class="ip-name">(gamma)</tspan></text>
    <text x="550" y="116" class="ip-desc">scales gene g's count in cell c, per class k</text>
    <text x="550" y="183" class="ip-glyph">&#952; <tspan class="ip-name">(theta)</tspan></text>
    <text x="550" y="203" class="ip-desc">scales cell c's total count, per class k</text>
    <text x="550" y="268" class="ip-glyph">&#951; <tspan class="ip-name">(eta)</tspan></text>
    <text x="550" y="288" class="ip-desc">scales gene g's count across all cells</text>
    <text x="550" y="353" class="ip-glyph" font-size="20">Inefficiency</text>
    <text x="550" y="373" class="ip-desc">scales the whole experiment at once</text>
  </g>
</svg>
<figcaption>The four factors by granularity.</figcaption>
</figure>

Four factors are needed because the mismatch occurs at all four levels at once: a
global scale difference between the cell type definitions and the experiment, a per-gene detection pattern, per-cell
variation in yield, and gene-by-cell noise. Each factor corrects the mismatch at its own
scale.

The factors are applied together. For each (cell, gene, class) triplet the expected
expression is the mean expression of the cell type definition times Inefficiency, eta,
theta and gamma.

The [demo](scale-factors-demo.md) puts three of them on sliders over five cells generated
from known class definitions.

## Inefficiencies

Each of eta, theta and gamma is the same statistic, a ratio of observed over expected:

$$
\text{factor} \;=\;
\frac{\text{observed}}{\text{expected}}
$$

- **gamma** compares the observed count of a *given gene in a given cell* with its
  expected count, under a candidate class.
- **theta** compares the observed total count of a *given cell* with its expected
  total, under a candidate class.
- **eta** compares the observed count of a *given gene across all cells* with its
  expected total.

A factor is greater than 1 when the observed count exceeds the expected count and less
than 1 otherwise. The factors differ only in the level of aggregation before the ratio
is formed.

## The priors

Each factor has a prior with mean 1, no rescaling, whose strength is a hyperparameter:
`rGene` for eta, `rTheta` for theta, `rSpot` for gamma. How far a factor moves from 1
depends on that hyperparameter against the number of counts available to estimate it,
and those counts differ by orders of magnitude between the three levels.

- **eta** is estimated from one gene's spots across the whole section, typically thousands. At
  the default `rGene` of 20 the prior is negligible and the data determine the posterior
  eta.
- **theta** is estimated from one cell's spots, typically tens. At the default `rTheta` of 25 the
  prior is comparable to the data. Lower it to around 2 and the data determine the posterior
  theta. It has to stay above 1: theta is
  $(\text{spots} + r_\theta - 1) / (r_\theta + \text{expected})$, so at 1 or below a cell with no
  spots gives a theta of zero or less. A low `rTheta` also lets a cell with very few spots be
  assigned to a real class on one or two spots;
  [the prior strength](../the-model/scale-factors.md#r-theta) describes the
  consequences.
- **gamma** is estimated from a given gene in a given cell, usually a fraction of a spot. The default
  `rSpot` of 2 therefore dominates, which is the intent: there is too little data at that level
  to estimate anything on its own. `rSpot` is also the dispersion of the negative binomial, since
  integrating gamma out is what produces it.

Raising a hyperparameter holds its factor near 1; lowering it lets the data set it. A
starting value for `rTheta` is the typical number of counts in a cell.

## The spatial factor (the MRF)

Cells of the same class tend to cluster in space, in layers or regions. When a cell is scored
against the classes, each class receives a bonus proportional to the sum of the
neighbours' probabilities of that class, closer neighbours weighted more. This is the
spatial term, a Markov random field (MRF), where a cell's label depends on its
neighbours' labels. `mrf_beta` sets its strength; at `0` only the gene counts and the
prior determine the class. It is described with the cell typing step in
[cell to class](cell-to-celltype.md#the-class-prior-and-the-spatial-term).

<!-- Parked, bring it back with the rest of the story:

A cell with almost no spots has nothing to weigh against its neighbours, so the spatial
term alone can decide its class.

On its own the sentence states the problem and leaves the reader with no way out, so it
needs the settings that hold such a cell back: a heavier prior on the Zero class through
`cell_type_weights`, and a larger `rTheta`, which keeps theta near 1 in a cell with few
counts instead of letting the data shrink it. Say when each one is the right knob, and
what it costs. -->

