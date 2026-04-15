# Class "QuadrupenFit"

Class "QuadrupenFit"

Class "QuadrupenFit"

## Details

Class of object returned by any fitting function of the quadrupen
package (`elastic.net` or `bounded.reg`).

This class comes with the usual
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), `show()`,
[`print()`](https://rdrr.io/r/base/print.html) and
[`deviance()`](https://rdrr.io/r/stats/deviance.html) S3 methods.

Specific R6 methods are available for model extraction
`QuadrupenFit$get_model()`, cross validation
`QuadrupenFit$cross_validate()`, stability selection
`QuadrupenFit$stability_path()`, criteria derivation
QuadrupenFit\$criteria() and plotting `QuadrupenFit$plot()`. They come
with equivalent S3 methods :
[`cross_validate()`](https://jchiquet.github.io/quadrupen/reference/cross_validate.md),
[`stability()`](https://jchiquet.github.io/quadrupen/reference/stability.md)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

The `"path"`plot is available as soon as a fit has been performed. For
the others, the appropriate post-treatments must have been made via the
methods `QuadrupenFit$criteria()`, `QuadrupenFit$cross_validate()` or
`QuadrupenFit$stability()`

All plots functions are given with the default arguments, except for
`labels` and `log_scale`. If you need more control, please use the
dedicated methods: `QuadrupenFit$plot_path()`,
[`InformationCriteria$plot()`](https://jchiquet.github.io/quadrupen/reference/InformationCriteria.md),
[`CrossValidation$plot()`](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md),
[`StabilityPath$plot()`](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
or the corresponding S3 methods.

Plot method for regularization path

## See also

See also
[`InformationCriteria`](https://jchiquet.github.io/quadrupen/reference/InformationCriteria.md),
[`CrossValidation`](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md)
and
[`StabilityPath`](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)

[`cross_validate()`](https://jchiquet.github.io/quadrupen/reference/cross_validate.md)
Stability selection for Quadrupen object

[`stability()`](https://jchiquet.github.io/quadrupen/reference/stability.md)
Penalized criteria based on estimation of degrees of freedom

[`criteria()`](https://jchiquet.github.io/quadrupen/reference/criteria.md)

## Active bindings

- `nvar`:

  number of coefficient (without intercept)

- `nobs`:

  sample size

- `dataModel`:

  an object with class
  [`DataModel`](https://jchiquet.github.io/quadrupen/reference/DataModel.md)
  storing the data

- `major_tuning`:

  vector of "leading" tuning parameters (either l1, linf or l2)

- `minor_tuning`:

  vector of "minor" tuning parameters (either l1 or l2)

- `optim_monitoring`:

  list monitoring the optimization

- `optim_config`:

  list with low level options used for optimization.

- `fitted`:

  Matrix of fitted values, each column corresponding to a value of
  `lambda1`.

- `coefficients`:

  Matrix (class `"dgCMatrix"`) of coefficients with respect to the
  original input. The number of rows corresponds the length of
  `lambda1`.

- `intercept`:

  A vector containing the successive values of the (unpenalized)
  intercept. Equals to zero if `intercept` has been set to `FALSE`.

- `debias`:

  logical, should we rely on the debias coefficient of the regularizer
  (if available) or not

- `residuals`:

  Matrix of residuals, each column corresponding to a value of
  `lambda1`.

- `deviance`:

  the model deviance

- `degrees_freedom`:

  Estimated degree of freedoms for the successive `lambda1`.

- `r_squared`:

  vector giving the coefficient of determination as a function of
  lambda1.

- `information_criteria`:

  object with class
  [`InformationCriteria`](https://jchiquet.github.io/quadrupen/reference/InformationCriteria.md)
  storing various information criteria (AIC, BIC, GCV, etc) for the
  current fit.

- `cross_validation`:

  object with class
  [`CrossValidation`](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md)
  storing output of CV job. Only available once method cross_validate
  has been called.

- `stability_path`:

  object with class
  [`StabilityPath`](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
  storing output of stability selection. Only available once method
  \$stability has been called.

## Methods

### Public methods

- [`QuadrupenFit$new()`](#method-QuadrupenFit-new)

- [`QuadrupenFit$show()`](#method-QuadrupenFit-show)

- [`QuadrupenFit$print()`](#method-QuadrupenFit-print)

- [`QuadrupenFit$fit()`](#method-QuadrupenFit-fit)

- [`QuadrupenFit$get_model()`](#method-QuadrupenFit-get_model)

- [`QuadrupenFit$predict()`](#method-QuadrupenFit-predict)

- [`QuadrupenFit$cross_validate()`](#method-QuadrupenFit-cross_validate)

- [`QuadrupenFit$stability()`](#method-QuadrupenFit-stability)

- [`QuadrupenFit$criteria()`](#method-QuadrupenFit-criteria)

- [`QuadrupenFit$plot()`](#method-QuadrupenFit-plot)

- [`QuadrupenFit$plot_path()`](#method-QuadrupenFit-plot_path)

- [`QuadrupenFit$clone()`](#method-QuadrupenFit-clone)

------------------------------------------------------------------------

### Method `new()`

Initialize a QuadrupenFit model

#### Usage

    QuadrupenFit$new(data, intercept, regParam)

#### Arguments

- `data`:

  a
  [DataModel](https://jchiquet.github.io/quadrupen/reference/DataModel.md)
  object

- `intercept`:

  a logical; should an intercept be included in the mode?

- `regParam`:

  a list with two elements, a vector and a scalar, for the
  regularization

------------------------------------------------------------------------

### Method `show()`

User friendly print method

#### Usage

    QuadrupenFit$show()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

User friendly print method

#### Usage

    QuadrupenFit$print()

------------------------------------------------------------------------

### Method `fit()`

function performing the optimization

#### Usage

    QuadrupenFit$fit(control)

#### Arguments

- `control`:

  list controlling the optimization process

------------------------------------------------------------------------

### Method `get_model()`

Model extraction

#### Usage

    QuadrupenFit$get_model(selection, type = c("coefficients", "penalty", "index"))

#### Arguments

- `selection`:

  either a character (model selection criteria) of a scalar (lambda
  value)

- `type`:

  character for the desired output

#### Returns

either a vector of coefficients, a scalar or the model index

------------------------------------------------------------------------

### Method [`predict()`](https://rdrr.io/r/stats/predict.html)

Predict response for new sample based on the current model

#### Usage

    QuadrupenFit$predict(newx = NULL, selection = NULL)

#### Arguments

- `newx`:

  matrix of new values for the regressor with which to predict. If
  omitted, the fitted values are used.

- `selection`:

  either a character (model selection criteria) of a scalar (lambda
  value)

#### Returns

a vector of predicted value Cross-validation for Quadrupen object

------------------------------------------------------------------------

### Method [`cross_validate()`](https://jchiquet.github.io/quadrupen/reference/cross_validate.md)

Function that computes K-fold cross-validated error of a `quadrupen`
fit, possibly on a grid of `lambda1`, `lambda2`.

#### Usage

    QuadrupenFit$cross_validate(
      K = 10,
      folds = split(sample(1:self$nobs), rep(1:K, length = self$nobs)),
      lambda2 = self$minor_tuning,
      verbose = TRUE,
      cores = 1
    )

#### Arguments

- `K`:

  integer indicating the number of folds. Default is 10.

- `folds`:

  list of `K` vectors that describes the folds to use for the
  cross-validation. By default, the folds are randomly sampled with the
  specified K. The same folds are used for each values of `lambda2`.

- `lambda2`:

  tunes the \\\ell_2\\-penalty (ridge-like) of the fit. If none is
  provided, a vector of values is generated and a CV is performed on a
  grid of `lambda2` and `lambda1`, using the same folds for each
  `lambda2`.

- `verbose`:

  logical; indicates if the progression (the current `lambda2` should be
  displayed. Default is `TRUE.`

- `cores`:

  the number of cores to use. The default uses 1 core (safer in case
  your BLAS/LAPACK libraries are multithreaded)

#### Returns

an object with class
[CrossValidation](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md)
is sent back and stored as a field of the original QuadrupenFit object.

------------------------------------------------------------------------

### Method [`stability()`](https://jchiquet.github.io/quadrupen/reference/stability.md)

Compute the stability path of a (possibly randomized) fitting procedure
as introduced by Meinshausen and Buhlmann (2010).

#### Usage

    QuadrupenFit$stability(
      n_subsamples = 50,
      subsample_size = floor(self$nobs/2),
      subsamples = replicate(n_subsamples, sample(1:self$nobs, subsample_size), simplify =
        FALSE),
      weakness = 1,
      verbose = TRUE,
      cores = 1
    )

#### Arguments

- `n_subsamples`:

  integer indicating the number of subsamplings used to estimate the
  selection probabilities. Default is 100.

- `subsample_size`:

  integer indicating the size of each subsamples. Default is
  `floor(n/2)`.

- `subsamples`:

  list with `subsamples` entries with vectors describing the folds to
  use for the stability procedure. By default, the folds are randomly
  sampled with the specified `n_subsamples` and `subsample_size`
  argument.

- `weakness`:

  Coefficient used for randomizing the weights of each features. Default
  is 1\` for no randomization. See details below.

- `verbose`:

  logical; indicates if the progression should be displayed. Default is
  `TRUE`.

- `cores`:

  the number of cores to use. The default uses 1 core (safer in case
  your BLAS/LAPACK libraries are multithreaded)

#### Returns

an object with class
[StabilityPath](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
is sent back and stored as a field of the original QuadrupenFit object.

------------------------------------------------------------------------

### Method [`criteria()`](https://jchiquet.github.io/quadrupen/reference/criteria.md)

Produce a plot or send back the values of some penalized criteria
accompanied with the vector(s) of parameters selected accordingly. The
default behavior plots the BIC and the AIC (with respective factor
\\\log(n)\\ and \\2\\) yet the user can specify any penalty.

#### Usage

    QuadrupenFit$criteria(
      penalty = setNames(c(2, log(self$nobs), log(self$nvar), log(self$nobs) + 2 *
        log(self$nvar)), c("AIC", "BIC", "mBIC", "eBIC")),
      sigma = NULL
    )

#### Arguments

- `penalty`:

  a vector with as many penalties a desired. The default contains the
  penalty corresponding to the AIC and the BIC (\\2\\ and \\\log(n)\\).
  Setting the "names" attribute, as done in the default definition,
  leads to outputs which are easier to read.

- `sigma`:

  scalar: an estimate of the residual variance. When available, it is
  plugged-in the criteria, which may be more relevant. If `NULL` (the
  default), it is estimated as usual (see details).

#### Returns

an object with class
[InformationCriteria](https://jchiquet.github.io/quadrupen/reference/InformationCriteria.md)
is sent back and stored as a field of the original QuadrupenFit object.

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot method for QuadrupenFit

#### Usage

    QuadrupenFit$plot(
      type = c("path", "criteria", "crossval", "stability"),
      log_scale = TRUE,
      labels = NULL
    )

#### Arguments

- `type`:

  the type of plot, either `"path"` for regularization path;
  `"criteria"` for BIC-like information criteria ; `"crossval"` for
  cross-validation plot ; and `"stability"` for stability path.

- `log_scale`:

  logical; indicates if a log-scale should be used when `xvar="lambda"`.
  Default is `TRUE`.

- `labels`:

  vector indicating the names associated to the plotted variables. When
  specified, a legend is drawn in order to identify each variable. Only
  relevant when the number of predictor is small. Remind that the
  intercept does not count. Default is `NULL`.

------------------------------------------------------------------------

### Method `plot_path()`

Produce a plot of the solution path of a QuadrupenFit object.

#### Usage

    QuadrupenFit$plot_path(
      xvar = c("lambda", "fraction", "df"),
      log_scale = TRUE,
      title = paste(self$penalty, " path", sep = ""),
      standardize = TRUE,
      labels = NULL
    )

#### Arguments

- `xvar`:

  variable to plot on the X-axis: either `"lambda"` (\\\ell_1\\ penalty
  level, or \\\ell_2\\ for ridge and \\\ell\_\infty\\) or `"fraction"`
  (\\\ell_1\\-norm of the coefficients) or `df` for estimated degrees of
  freedom. Default is set to `"lambda"`.

- `log_scale`:

  logical; indicates if a log-scale should be used when `xvar="lambda"`.
  Default is `TRUE`.

- `title`:

  the title. Default is set to the model name followed by what is on the
  Y-axis.

- `standardize`:

  logical; standardize the coefficients before plotting (with the norm
  of the predictor). Default is `TRUE`.

- `labels`:

  vector indicating the names associated to the plotted variables. When
  specified, a legend is drawn in order to identify each variable. Only
  relevant when the number of predictor is small. Remind that the
  intercept does not count. Default is `NULL`.

#### Returns

a ggplot2 object .

#### Examples

    \dontrun{
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

    ## Plot the Lasso path
    plot(lasso(x,y), title="Lasso solution path")
    ## Plot the Elastic-net path
    plot(elastic.net(x,y), title = "Elastic-net solution path")
    ## Plot the Elastic-net path (fraction on X-axis, unstandardized coefficient)
    plot(elastic.net(x,y, lambda2=10), standardize=FALSE, xvar="fraction")
    ## Plot the Bounded regression path (fraction on X-axis)
    plot(bounded.reg(x,y, lambda2=10), xvar="fraction")
    }

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    QuadrupenFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
## ------------------------------------------------
## Method `QuadrupenFit$plot_path`
## ------------------------------------------------

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

## Plot the Lasso path
plot(lasso(x,y), title="Lasso solution path")
## Plot the Elastic-net path
plot(elastic.net(x,y), title = "Elastic-net solution path")
## Plot the Elastic-net path (fraction on X-axis, unstandardized coefficient)
plot(elastic.net(x,y, lambda2=10), standardize=FALSE, xvar="fraction")
## Plot the Bounded regression path (fraction on X-axis)
plot(bounded.reg(x,y, lambda2=10), xvar="fraction")
} # }
```
