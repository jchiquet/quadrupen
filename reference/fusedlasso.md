# A function for fitting generalized fused-Lasso problems

This function fits the standard version of the fused lasso. It can take
a general matrix x and provides for possible weights on the `lambda1`
and `lambda2` penalties.

## Usage

``` r
fusedlasso(
  x,
  y,
  lambda1 = NULL,
  lambda2 = 1,
  pen_fused = c("L1", "L2", "Huber"),
  penscale = rep(1, ncol(x)),
  struct = NULL,
  intercept = TRUE,
  normalize = TRUE,
  nlambda1 = ifelse(is.null(lambda1), 50, length(lambda1)),
  minratio = 0.01,
  maxfeat = ifelse(lambda2 < 1, min(nrow(x), ncol(x)), min(2 * nrow(x), ncol(x))),
  beta0 = rep(0, ncol(x)),
  control = list()
)
```

## Arguments

- x:

  matrix of features, possibly sparsely encoded (experimental). Do NOT
  include intercept. When normalized os `TRUE`, coefficients will then
  be rescaled to the original scale.

- y:

  response vector.

- lambda1:

  sequence of decreasing \\\ell_1\\-penalty levels. If `NULL` (the
  default), a vector is generated with `nlambda1` entries, starting from
  a guessed level `lambda1.max` where only the intercept is included,
  then shrunken to `minratio*lambda1.max`.

- lambda2:

  real scalar; tunes the \\\ell_2\\ penalty in the Elastic-net. Default
  is 0.01. Set to 0 to recover the Lasso.

- pen_fused:

  penalty used for fusing the variables (either L1, L2 or Huber).
  Default is L1

- penscale:

  vector with real positive values that weight the \\\ell_1\\-penalty of
  each feature. Default set all weights to 1.

- struct:

  Description of the graph that corresponds to the `lambda2` penalty
  structure. If `NULL` (the default) a chain graph is assumed, like in
  the standard fused-lasso. If a matrix is given, interpreted as a
  symmetric adjacency matrix

- intercept:

  logical; indicates if an intercept should be included in the model.
  Default is `TRUE`.

- normalize:

  logical; indicates if variables should be normalized to have unit L2
  norm before fitting. Default is `TRUE`.

- nlambda1:

  integer that indicates the number of values to put in the `lambda1`
  vector. Ignored if `lambda1` is provided.

- minratio:

  minimal value of \\\ell_1\\-part of the penalty that will be tried, as
  a fraction of the maximal `lambda1` value. A too small value might
  lead to unstability at the end of the solution path corresponding to
  small `lambda1` combined with \\\lambda_2=0\\. The default value tries
  to avoid this, adapting to the '\\n\<p\\' context. Ignored if
  `lambda1` is provided.

- maxfeat:

  integer; limits the number of features ever to enter the model; i.e.,
  non-zero coefficients for the Elastic-net: the algorithm stops if this
  number is exceeded and `lambda1` is cut at the corresponding level.
  Default is `min(nrow(x),ncol(x))` for small `lambda2` (\<0.01) and
  `min(4*nrow(x),ncol(x))` otherwise. Use with care, as it considerably
  changes the computation time.

- beta0:

  a starting point for the vector of parameter. By default, will
  initialized zero. May save time in some situation.

- control:

  list of argument controlling low level options of the algorithm:

  - `verbose`: logical; ener verbose mode

  - `timer`: logical; use to record the timing of the algorithm. Default
    is `FALSE`.

  - `maxiterin` Maximum number of iterations in the inner loop to run.

  - `maxiterout` Maximum number of iterations in the outer loop to run.

  - `maxactivation` Maximum number of previously inactive variables to
    activate at the same time

  - `accuracy` Accuracy at which the algorithm will stop.

  - `fusioncheck` Should the fused sets be checked for breaking up?

  - `verbose` Should the function give some output what it is doing?

## Value

an object with class
[FusedLassoFit](https://jchiquet.github.io/quadrupen/reference/FusedLassoFit.md),
inheriting from
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md).

## Details

The optimized criterion is the following:

β^(hat)_(λ₁,λ₂) = argmin_(β) 1/2 RSS(&beta) + λ₁ \| D β \|₁ + λ/2 ₂
β^(T) S β,

where \\D\\ is a diagonal matrix, whose diagonal terms are provided as a
vector by the `penscale` argument. The \\\ell_1\\ fusion penalty is
structured by a possibily weighted graph \\G\\ provided via the `struct`
argument, as a symmetrix (undirected) adjacency matrix.

## Author

Original code by Holger Hoefling, refactoring by Julien Chiquet

## Examples

``` r
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor <- 0.75
Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- Matrix::bdiag(Soo,Sww,Soo,Sww,Soo)
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

res <- fusedlasso(x, y, lambda2=5)
G <- igraph::make_ring(ncol(x)) |> igraph::as_adjacency_matrix(sparse = FALSE)
resG <- fusedlasso(x, y, lambda2=5, struct = G)
plot(res)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'drop': error in evaluating the argument 'x' in selecting a method for function 'colSums': object 'x' not found
plot(resG)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'drop': error in evaluating the argument 'x' in selecting a method for function 'colSums': object 'x' not found
```
