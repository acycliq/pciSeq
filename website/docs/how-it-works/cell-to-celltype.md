
# 3. Assigning cells to classes

Every cell is scored against every class by how well its gene counts agree with
the class's warped expected expression, and the scores are converted to probabilities.
The inputs are the gene counts per cell from the current spot assignment and the warped
definitions from the previous step.

## The expression score

The score of a cell for a class is the log probability of the cell's counts under that
class's expected expression, summed over genes. A class that predicts genes the cell
lacks, or omits genes the cell has, scores low. A few strong marker genes can dominate
the sum, which is what separates closely related classes.

## The softmax

A softmax over the scores gives each cell a distribution over the classes, for example
0.80 on one class, 0.15 on a close relative and 0.05 elsewhere. A confident assignment
concentrates the mass on one class, an ambiguous one spreads it.

## The class prior and the spatial term

Two further terms are added to the score.

- **The class prior.** With `cell_type_prior` set to `uniform`, the default, the weights
  are those of `cell_type_weights` and stay fixed; by default every real class has the same
  weight. With `weighted` the weights are estimated from the
  data with a [Dirichlet update](../the-model/class-prior.md), so common classes are
  favoured over rare ones.

- **The spatial term (the MRF).** Cells of the same class cluster in space. Each class
  receives a bonus proportional to the weighted number of the cell's neighbours that
  carry it, with closer neighbours weighted more, on a scale set per cell from its own
  neighbour distances. The term is added to the score, so strong gene evidence can
  outweigh it, and it is bounded. Its strength is `mrf_beta`; `0` switches it off. The
  weights are derived in
  [the model](../the-model/cell-class.md#weighting-the-neighbours-by-distance).
  Sister classes can be pooled so the neighbours support all of them equally and the
  expression alone decides between them, see
  [pooling sister classes](../the-model/cell-class.md#pooling-sister-classes).

## The Zero class

One class, Zero, expects no expression. It takes cells with zero or few reads that no
real class explains; such a cell can be a real cell of a type the gene panel does not
target.
Without it such cells would be forced onto a real class. A near-empty cell has no counts
to weigh against its neighbours, so the spatial term alone can determine its
class.

The output, a distribution over classes for every cell, is the input to
[spot assignment](spots-to-cells.md).
