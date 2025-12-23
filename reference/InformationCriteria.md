# Class InformationCriteria

Class InformationCriteria

Class InformationCriteria

## Details

Class of object returned by the
[`QuadrupenFit$criteria()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
method or the
[`criteria()`](https://jchiquet.github.io/quadrupen/reference/criteria.md)
function. Owns [`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## Active bindings

- `data`:

  a data frame containing the values of various information criteria
  (AIC, BIC, lmBIC, eBIC, GCV) along the value of `lambda1`.

- `lambda`:

  vector of \\\lambda_1\\ (\\\ell_1\\ or \\\ell\_\infty\\ penalty
  levels) for which each cross-validation has been performed.

- `names`:

  a vector of characters storing the names of the precomputed criteria

## Methods

### Public methods

- [`InformationCriteria$new()`](#method-InformationCriteria-new)

- [`InformationCriteria$show()`](#method-InformationCriteria-show)

- [`InformationCriteria$print()`](#method-InformationCriteria-print)

- [`InformationCriteria$plot()`](#method-InformationCriteria-plot)

- [`InformationCriteria$clone()`](#method-InformationCriteria-clone)

------------------------------------------------------------------------

### Method `new()`

Constructor for a InformationCriteria object Should be called internally
by an object
[`QuadrupenFit$criteria()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

#### Usage

    InformationCriteria$new(value)

#### Arguments

- `value`:

  data frame storing output of
  [`QuadrupenFit$criteria()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

------------------------------------------------------------------------

### Method `show()`

User friendly print method

#### Usage

    InformationCriteria$show()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

User friendly print method

#### Usage

    InformationCriteria$print()

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Plot the the desired criteria

#### Usage

    InformationCriteria$plot(
      criteria = self$names,
      log_scale = TRUE,
      xvar = c("lambda", "fraction", "df"),
      title = "Information Criteria"
    )

#### Arguments

- `criteria`:

  a vector of character with the criteria to plot. The default plot all
  the criteria available (stored in the field `names`)

- `log_scale`:

  logical; indicates if a log-scale should be used when `xvar="lambda"`.
  Default is `TRUE`.

- `xvar`:

  variable to plot on the X-axis: either `"df"` (the estimated degrees
  of freedom), `"lambda"` (\\\lambda_1\\ penalty level) or `"fraction"`
  (\\\ell_1\\-norm of the coefficients). Default is set to `"lambda"`.

- `title`:

  graph title

#### Returns

a ggplot2 object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    InformationCriteria$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
