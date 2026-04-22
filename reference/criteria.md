# Penalized criteria based on estimation of degrees of freedom

Produce a plot or send back the values of some penalized criteria
accompanied with the vector(s) of parameters selected accordingly. The
default behavior plots the BIC and the AIC (with respective factor
\\\log(n)\\ and \\2\\) yet the user can specify any penalty.

## Usage

``` r
criteria(
  object,
  penalty = setNames(c(2, log(object$nobs), log(object$nvar), log(object$nobs) + 2 *
    log(object$nvar)), c("AIC", "BIC", "mBIC", "eBIC")),
  sigma = NULL
)

# S3 method for class 'QuadrupenFit'
criteria(
  object,
  penalty = setNames(c(2, log(object$nobs), log(object$nvar), log(object$nobs) + 2 *
    log(object$nvar)), c("AIC", "BIC", "mBIC", "eBIC")),
  sigma = NULL
)
```

## Arguments

- object:

  output of a fitting procedure of the quadrupen package (e.g.
  [`elastic.net()`](https://jchiquet.github.io/quadrupen/reference/elastic.net.md)).

- penalty:

  a vector with as many penalties a desired. The default contains the
  penalty corresponding to the AIC and the BIC (\\2\\ and \\\log(n)\\).
  Setting the "names" attribute, as done in the default definition,
  leads to outputs which are easier to read.

- sigma:

  scalar: an estimate of the residual variance. When available, it is
  plugged-in the criteria, which may be more relevant. If `NULL` (the
  default), it is estimated as usual (see details).

## Value

an object with class
[InformationCriteria](https://jchiquet.github.io/quadrupen/reference/InformationCriteria.md)
is sent back and stored as a field of the original
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
object.

## Methods (by class)

- `criteria(QuadrupenFit)`: S3 method for information criteria of a
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

## Note

When `sigma` is provided, the criterion takes the form

RSS + penalty \* df / n \* sigma²

When it is unknown, it writes

n\*log(RSS) + penalty \* df

Estimation of the degrees of freedom (for the elastic-net, the LASSO and
also bounded regression) are computed by applying and adapting the
results of Tibshirani and Taylor (see references below).

## References

Ryan Tibshirani and Jonathan Taylor. Degrees of freedom in lasso
problems, Annals of Statistics, 40(2) 2012.

## Examples

``` r
if (FALSE) { # \dontrun{
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor <- 0.75
Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

## Plot penalized criteria for the Elastic-net path
criteria(elastic.net(x,y, lambda2=1))

#' Plot penalized criteria for the Bounded regression
criteria(bounded.reg(x,y, lambda2=1))
} # }
```
