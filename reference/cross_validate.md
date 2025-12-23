# Cross-validation for Quadrupen object

Function that computes K-fold cross-validated error of a `quadrupen`
fit, possibly on a grid of `lambda1`, `lambda2`.

## Usage

``` r
cross_validate(
  object,
  K = 10,
  folds = split(sample(1:object$nobs), rep(1:K, length = object$nobs)),
  lambda2 = object$minor_tuning,
  verbose = TRUE,
  cores = parallel::detectCores() - 2
)

# S3 method for class 'QuadrupenFit'
cross_validate(
  object,
  K = 10,
  folds = split(sample(1:object$nobs), rep(1:K, length = object$nobs)),
  lambda2 = object$minor_tuning,
  verbose = TRUE,
  cores = parallel::detectCores() - 2
)
```

## Arguments

- object:

  an R6 object with class
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

- K:

  integer indicating the number of folds. Default is 10.

- folds:

  list of `K` vectors that describes the folds to use for the
  cross-validation. By default, the folds are randomly sampled with the
  specified K. The same folds are used for each values of `lambda2`.

- lambda2:

  tunes the \\\ell_2\\-penalty (ridge-like) of the fit. If none is
  provided, a vector of values is generated and a CV is performed on a
  grid of `lambda2` and `lambda1`, using the same folds for each
  `lambda2`.

- verbose:

  logical; indicates if the progression (the current `lambda2` should be
  displayed. Default is `TRUE.`

- cores:

  the number of cores to use. The default uses all the cores available.

## Value

an object with class
[`CrossValidation`](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md)
is sent back and stored as a field of the original
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
object.

## Methods (by class)

- `cross_validate(QuadrupenFit)`: S3 method for cross-validation of a
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

## Note

If the user runs the fitting method with option `'bulletproof'` set to
`FALSE`, the algorithm may stop at an early stage of the path. Early
stops are handled internally, in order to provide results on the same
grid of penalty tuned by \\\lambda_1\\. This is done by means of `NA`
values, so as mean and standard error are consistently evaluated. If,
while cross-validating, the procedure experiences too much early
stoppings, a warning is sent to the user, in which case you should
reconsider the grid of `lambda1` used for the cross-validation. If
`bulletproof` is `TRUE` (the default), there is nothing to worry about,
except a possible slow down when any switching to the proximal algorithm
is required.

## Examples

``` r
if (FALSE) { # \dontrun{
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor  <- 0.75
Soo  <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variable
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.1
diag(Sigma) <- 1
n <- 100
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

enet <- elastic.net(x, y, nlambda1=50)

## Use fewer lambda1 values by overwritting the default parameters
## and cross-validate over the sequences lambda1 and lambda2
cv.grid <- cross_validate(enet, lambda2=10^seq(2,-2,len=50))
## Rerun simple cross-validation with the appropriate lambda2
cv.10K <- crossval(x,y, lambda2=cv.grid$lambda2_min)
## Try leave one out also
cv.loo <- crossval(x,y, K=n, lambda2=cv.grid$lambda2_min)

plot(cv.grid)
plot(cv.10K)
plot(cv.loo)

## Performance for selection purpose
cat("\nFalse positives with the minimal 10-CV choice: ", sum(sign(beta) != sign(cv.10K$beta_min )))
cat("\nFalse positives with the minimal LOO-CV choice: ", sum(sign(beta) != sign(cv.loo$beta_min)))
} # }
```
