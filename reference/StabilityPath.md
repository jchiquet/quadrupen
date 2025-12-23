# Class StabilityPath

Class StabilityPath

Class StabilityPath

## Details

Class of object returned by the
[`QuadrupenFit$cross_validate()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
method or the
[`cross_validate()`](https://jchiquet.github.io/quadrupen/reference/cross_validate.md)
function. Owns [`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## Public fields

- `probabilities`:

  a `Matrix` object containing the estimated probabilities of selection
  along the path of solutions.

- `regParam`:

  a lsit with the levels of the regularizing parameters used

- `subsamples`:

  a list that contains the folds used for each subsample.

## Active bindings

- `nvar`:

  number of variables (without intercept)

- `nobs`:

  number of observation/sample size

- `nonzero`:

  variables with a non-null probability of selection along the stability
  path

- `nonzeroprob`:

  subset of the probabilityes stability path on the nonzero variables

## Methods

### Public methods

- [`StabilityPath$new()`](#method-StabilityPath-new)

- [`StabilityPath$show()`](#method-StabilityPath-show)

- [`StabilityPath$print()`](#method-StabilityPath-print)

- [`StabilityPath$selection()`](#method-StabilityPath-selection)

- [`StabilityPath$plot()`](#method-StabilityPath-plot)

- [`StabilityPath$clone()`](#method-StabilityPath-clone)

------------------------------------------------------------------------

### Method `new()`

Constructor for a StabilityPath object Should be called internally by an
object
[`QuadrupenFit$stability()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

#### Usage

    StabilityPath$new(probabilities, regParam, subsamples)

#### Arguments

- `probabilities`:

  a `Matrix` object containing the estimated probabilities of selection
  along the path of solutions.

- `regParam`:

  a lsit with the levels of the regularizing parameters used

- `subsamples`:

  a list that contains the folds used for each subsample.

------------------------------------------------------------------------

### Method `show()`

User friendly print method

#### Usage

    StabilityPath$show()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

User friendly print method

#### Usage

    StabilityPath$print()

------------------------------------------------------------------------

### Method `selection()`

Perform variable selection based on the stability path

#### Usage

    StabilityPath$selection(
      sel_mode = c("rank", "PFER"),
      cutoff = 0.75,
      PFER = 2,
      nvarsel = floor(self$nobs/log(self$nvar))
    )

#### Arguments

- `sel_mode`:

  a character string, either `"rank"` or `"PFER"`. In the first case,
  the selection is based on the rank of total probabilties by variables
  along the path: the first `nvarsel` variables are selected (see
  below). In the second case, the PFER control is used as described in
  Meinshausen and Buhlmannn's paper. Default is `"rank"`.

- `cutoff`:

  value of the cutoff probability (only relevant when `sel_mode` equals
  `"PFER"`).

- `PFER`:

  value of the per-family error rate to control (only relevant when
  `sel_mode` equals `"PFER"`).

- `nvarsel`:

  number of variables selected (only relevant when `sel_mode` equals
  `"rank"`. Default is `floor(n/log(p))`.

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Produce a plot of the stability path obtained by stability selection.

#### Usage

    StabilityPath$plot(
      xvar = "lambda",
      title = "Stability path",
      labels = rep("unknown status", p),
      sel_mode = c("rank", "PFER"),
      cutoff = 0.75,
      PFER = 2,
      nvarsel = floor(self$nobs/log(self$nvar))
    )

#### Arguments

- `xvar`:

  variable to plot on the X-axis: either `"lambda"` (first penalty
  level) or `"fraction"` (fraction of the penalty level applied tune by
  \\\lambda_1\\). Default is `"lambda"`.

- `title`:

  title title. If none given, a somewhat appropriate title is
  automatically generated.

- `labels`:

  an optional vector of labels for each variable in the path (e.g.,
  'relevant'/'irrelevant'). See examples.

- `sel_mode`:

  a character string, either `"rank"` or `"PFER"`. In the first case,
  the selection is based on the rank of total probabilties by variables
  along the path: the first `nvarsel` variables are selected (see
  below). In the second case, the PFER control is used as described in
  Meinshausen and Buhlmannn's paper. Default is `"rank"`.

- `cutoff`:

  value of the cutoff probability (only relevant when `sel_mode` equals
  `"PFER"`).

- `PFER`:

  value of the per-family error rate to control (only relevant when
  `sel_mode` equals `"PFER"`).

- `nvarsel`:

  number of variables selected (only relevant when `sel_mode` equals
  `"rank"`. Default is `floor(n/log(p))`.

- `log_scale`:

  logical; indicates if a log-scale should be used when `xvar="lambda"`.
  Default is `TRUE`.

- `plot`:

  logical; indicates if the graph should be plotted. Default is `TRUE`.
  If `FALSE`, only the ggplot2 object is sent back.

#### Returns

a list with a ggplot2 object which can be plotted via the `print`
method, and a vector of selected variables corresponding to method of
choice (`"rank"` or `"PFER"`).

#### Examples

    \dontrun{
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

    ## Call to stability selection function, 200 subsampling
    stab <- stability(x,y, subsamples=200, lambda2=1, minratio=1e-2)

    ## Build the plot an recover the selected variable
    plot(stab, labels=labels)
    plot(stab, xvar="fraction", labels=labels, sel_mode="PFER", cutoff=0.75, PFER=2)
    }

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    StabilityPath$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
## ------------------------------------------------
## Method `StabilityPath$plot`
## ------------------------------------------------

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

## Call to stability selection function, 200 subsampling
stab <- stability(x,y, subsamples=200, lambda2=1, minratio=1e-2)

## Build the plot an recover the selected variable
plot(stab, labels=labels)
plot(stab, xvar="fraction", labels=labels, sel_mode="PFER", cutoff=0.75, PFER=2)
} # }
```
