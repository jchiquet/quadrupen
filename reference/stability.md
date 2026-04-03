# Stability selection for Quadrupen object

Compute the stability path of a (possibly randomized) fitting procedure
as introduced by Meinshausen and Buhlmann (2010).

## Usage

``` r
stability(
  object,
  n_subsamples = 50,
  subsample_size = floor(object$nobs/2),
  subsamples = replicate(n_subsamples, sample(1:object$nobs, subsample_size), simplify =
    FALSE),
  weakness = 1,
  verbose = TRUE,
  cores = parallel::detectCores() - 2
)

# S3 method for class 'QuadrupenFit'
stability(
  object,
  n_subsamples = 50,
  subsample_size = floor(object$nobs/2),
  subsamples = replicate(n_subsamples, sample(1:object$nobs, subsample_size), simplify =
    FALSE),
  weakness = 1,
  verbose = TRUE,
  cores = parallel::detectCores() - 2
)
```

## Arguments

- object:

  an R6 object with class
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

- n_subsamples:

  integer indicating the number of subsamplings used to estimate the
  selection probabilities. Default is 100.

- subsample_size:

  integer indicating the size of each subsamples. Default is
  `floor(n/2)`.

- subsamples:

  list with `subsamples` entries with vectors describing the folds to
  use for the stability procedure. By default, the folds are randomly
  sampled with the specified `n_subsamples` and `subsample_size`
  argument.

- weakness:

  Coefficient used for randomizing the weights of each features. Default
  is 1\` for no randomization. See details below.

- verbose:

  logical; indicates if the progression should be displayed. Default is
  `TRUE`.

- cores:

  the number of cores to use. The default uses all the cores available.

## Value

an object with class
[`StabilityPath`](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
is sent back and stored as a field of the original
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
object.

## Methods (by class)

- `stability(QuadrupenFit)`: S3 method for stability selection of a
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

## Note

When `weakness < 1`, the `penscale` argument that weights the penalty
tuned by \\\lambda_1\\ is perturbed (divided) for each subsample by a
random variable uniformly distributed on \[&#945;,1\], where &#945; is
the weakness parameter.

If the user runs the fitting method with option `'bulletproof'` set to
`FALSE`, the algorithm may stop at an early stage of the path. Early
stops of the underlying fitting function are handled internally, in the
following way: we chose to simply skip the results associated with such
runs, in order not to bias the stability selection procedure. If it
occurs too often, a warning is sent to the user, in which case you
should reconsider the grid of `lambda1` for stability selection. If
`bulletproof` is `TRUE` (the default), there is nothing to worry about,
except a possible slow down when any switching to the proximal algorithm
is required.

## References

N. Meinshausen and P. Buhlmann (2010). Stability Selection, JRSS(B).

## Examples

``` r
if (FALSE) { # \dontrun{
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
diag(Sigma) <- 1
n <- 100
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

## Build a vector of label for true nonzeros
labels <- rep("irrelevant", length(beta))
labels[beta != 0] <- c("relevant")
labels <- factor(labels, ordered=TRUE, levels=c("relevant","irrelevant"))

enet <- elastic.net(x, y, lambda2 = 10, struct = solve(Sigma), minratio = 1e-2)
stab <- stability(enet, n_subsamples = 200)

## Build the plot an recover the selected variable
plot(stab, labels=labels)
stabpath <- plot(stab, xvar="fraction", labels=labels, sel_mode="PFER", cutoff=0.75, PFER=1)

cat("\nFalse positives for the randomized Elastic-net with stability selection: ",
     sum(labels[stab$selection()] != "relevant"))
cat("\nDONE.\n")
} # }
```
