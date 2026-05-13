# Sparsity by Worst-Case Quadratic Penalties

Fits the solution paths of classical sparse regression models with
efficient active set algorithms by solving small sub-problems. Depending
on the penalty, the sub-problems can be solved exactly (*i.e.* for the
LASSO) or with generic solvers. The available optimizer includes
quadratic solvers, Newton-based approaches and generic FISTA or PGD
algorithms. Also provides a few methods for model selection purpose
(information criteria, cross-validation, stability selection).

## Details

**Quadrupen** covers the following regularizers

- LASSO (Least Absolute Shrinkage and Selection Operator)

- SCAD (Smoothly Clip Absolute Deviation)

- MCP (Minimax Concave Penalty)

- Group-LASSO (L1/L2 or L1/Linfty)

- Cooperative-LASSO

- Sparse Group-LASSO and Sparse Cooperative-LASSO

- Bounded Regression (L-infty norm).

For all these regularizers, **Quadrupen** offers the possibility to add
an ridge-like "structured" penalty to embed some external knowledge
about the statistical dependence between the features. This is sometimes
referred to as the "Structured Elastic-Net".

We also provide in the package the implementation of the Generalized
Fused-LASSO originally proposed by Holger Höfling now archived from CRAN
([original repo here](https://github.com/cran/FusedLasso)).

While likely not as fast as highly specialized packages like *glmnet*,
the use of a working set algorithm combined with efficient solvers,
sparse matrix support when applicable, and templated C++ code makes it
both competitive and versatile.

### Features:

The more important functions of the package are
[`sparse_lm()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md)
and
[`group_sparse_lm()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)
functions, which fits a linear model either with element-wise or
group-wise sparsity. The functions
[`lasso()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md),
[`elastic_net()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md),
[`scad()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md),
[`mcp()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md)
and
[`group_lasso()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md),
[`group_l1linf()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md),
[`coop_lasso()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md),
[`sparse_group_lasso()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md),
[`sparse_coop_lasso()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)
are only aliases for these two main functions.

The functions
[`lava()`](https://jchiquet.github.io/quadrupen/reference/lava.md),
[`group_lava()`](https://jchiquet.github.io/quadrupen/reference/group_lava.md),
[`fused_lasso()`](https://jchiquet.github.io/quadrupen/reference/fused_lasso.md)
and
[`bounded_reg()`](https://jchiquet.github.io/quadrupen/reference/bounded_reg.md)
are also available and specific.

We also included R6 and S3 methods for plotting, cross-validation and
for the stability selection procedure of Meinshausen and Buhlmann
(2010).

### Algorithm:

The general strategy of the algorithm relies on maintaining an active
set of variables, starting from a vector of zeros. The underlying
optimization problem is solved only on the activated variables, thus
handling with small smooth problems with increasing size. Hence, by
considering a decreasing grid of values for the penalty \\\lambda_1\\
and fixing \\\lambda_2\\, we may explore the whole path of solutions at
a reasonable numerical cost, providing that \\\lambda_1\\ does not end
up too small.

For the \\\ell_1\\-based methods (available in the `elastic_net`
function), the size of the underlying problems solved is related to the
number of nonzero coefficients in the vector of parameters. With the
\\\ell\_\infty\\-norm, (available in the `boundary.reg` function), we do
not produce sparse estimator. Nevertheless, the size of the systems
solved along the path deals with the number of unbounded variables for
the current penalty level, which is quite smaller than the number of
predictors for a reasonable \\\lambda_1\\. The same kind of proposal was
made in Zhao, Rocha and Yu (2009).

Underlying optimization is performed by direct resolution of (quadratic)
sub problems, which is the main purpose of this package. We also
implemented the popular and versatile proximal (FISTA) approaches for
routine checks and numerical comparisons. A Proximal Gradient Descent
approach with Anderson acceleration is also included.

The default setting uses the most appropriate solver (quadratic or
FISTA). The quadratic approach, which gave its name to the package, has
been optimized to be the method of choice for small and medium scale
problems, and produce very accurate solutions (in particular for
Elastic-Net/Lasso, Group-Lasso and Bounded Regression). However, the
first order methods (PGD and FISTA) remain competitive in particular in
situations where the problem is close to singular, in which case the
Cholesky/Eigen value decomposition used in the quadratic solver can be
computationally unstable.

## References

- Yves Grandvalet, Julien Chiquet and Christophe Ambroise, [Sparsity by
  Worst-case Quadratic Penalties](http://arxiv.org/abs/1210.2077), arXiv
  preprint, 2012.

- Fan, Jianqing, and Runze Li. “Variable selection via nonconcave
  penalized likelihood and its oracle properties.” JASA, 2001

- Nicolas Meinshausen and Peter Buhlmann. Stability Selection, JRSS(B),
  2010.

- Hoefling, Holger. “A path algorithm for the fused lasso signal
  approximator.” JCGS, 2010.

- Martin Slawski, Wolfgang zu Castell, and Gerhard Tutz. Feature
  selection guided by structural information, AOAS,

1.  

- Zhang, C. H. Nearly unbiased variable selection under minimax concave
  penalty. The Annals of Statistics, 2010.

- Yuan, Ming, and Yi Lin. “Model selection and estimation in regression
  with grouped variables.”, JRSS(B), 2006.

- Simon, Noah, et al. “A sparse-group lasso.” JCGS, 2013

- Peng Zhao, Guillerme Rocha and Bin Yu. The composite absolute
  penalties family for grouped and hierarchical variable selection, The
  Annals of Statistics, 2009.

- Hui Zou and Trevor Hastie. Regularization and variable selection via
  the elastic net, JRSS(B), 2006.

- Robert Tibshirani. Regression Shrinkage and Selection via the Lasso,
  JRSS(B), 1996.

## See also

Useful links:

- <https://github.com/jchiquet/quadrupen>

- <https://jchiquet.github.io/quadrupen/>

- Report bugs at <https://github.com/jchiquet/quadrupen/issues>

## Author

Julien Chiquet <julien.chiquet@inrae.fr>
