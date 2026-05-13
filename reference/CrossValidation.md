# Class CrossValidation

Class of object returned by the
[`QuadrupenFit$cross_validate()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
method or the
[`cross_validate()`](https://jchiquet.github.io/quadrupen/reference/cross_validate.md)
function. Owns [`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## Active bindings

- `data`:

  a data frame containing the mean cross-validated error and its
  associated standard error for each values of `lambda1` and `lambda2`.

- `folds`:

  list of `K` vectors indicating the folds used for cross-validation.

- `lambda1`:

  vector of \\\lambda_1\\ (\\\ell_1\\ or \\\ell\_\infty\\ penalty
  levels) for which each cross-validation has been performed.

- `lambda2`:

  vector (or scalar) of \\\ell_2\\-penalty levels for which each
  cross-validation has been performed.

- `lambda1_min`:

  level of \\\lambda_1\\ that minimizes the error estimated by
  cross-validation.

- `lambda2_min`:

  level of \\\lambda_2\\ that minimizes the error estimated by
  cross-validation.

- `lambda1_1se`:

  largest level of \\\lambda_1\\ such as the cross-validated error is
  within 1 standard error of the minimum.

- `lambda2_1se`:

  largest level of \\\lambda_2\\ such that the cross-validated error is
  within 1 standard error of the minimum (only relevant for ridge
  regression).

## Methods

### Public methods

- [`CrossValidation$new()`](#method-CrossValidation-initialize)

- [`CrossValidation$show()`](#method-CrossValidation-show)

- [`CrossValidation$print()`](#method-CrossValidation-print)

- [`CrossValidation$plotCV_1D()`](#method-CrossValidation-plotCV_1D)

- [`CrossValidation$plotCV_2D()`](#method-CrossValidation-plotCV_2D)

- [`CrossValidation$plot()`](#method-CrossValidation-plot)

- [`CrossValidation$clone()`](#method-CrossValidation-clone)

------------------------------------------------------------------------

### `CrossValidation$new()`

Constructor for a CrossValidation object Should be called internally by
an object
[`QuadrupenFit$cross_validate()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

#### Usage

    CrossValidation$new(cv_error, folds)

#### Arguments

- `cv_error`:

  data frame storing output of a cv job

- `folds`:

  list of K folds used for cross-validation

------------------------------------------------------------------------

### `CrossValidation$show()`

User friendly print method

#### Usage

    CrossValidation$show()

------------------------------------------------------------------------

### `CrossValidation$print()`

User friendly print method

#### Usage

    CrossValidation$print()

------------------------------------------------------------------------

### `CrossValidation$plotCV_1D()`

Plot 1-dimensional cross-validation

#### Usage

    CrossValidation$plotCV_1D(log_scale = TRUE, title = "Cross-validation error")

#### Arguments

- `log_scale`:

  logical, should a log-scale be used for the x-axis

- `title`:

  graph title

#### Returns

a ggplot object

------------------------------------------------------------------------

### `CrossValidation$plotCV_2D()`

Plot 2-dimensional cross-validation output (grid lambda1 x lambda2)

#### Usage

    CrossValidation$plotCV_2D(title = "Cross-validation error")

#### Arguments

- `title`:

  graph title

#### Returns

a ggplot2 object

------------------------------------------------------------------------

### `CrossValidation$plot()`

Plot cross-validation job by choosing the most appropriate output (1D-
or 2D)

#### Usage

    CrossValidation$plot(log_scale = TRUE, title = "Cross-validation error")

#### Arguments

- `log_scale`:

  logical, should a log-scale be used for the x-axis

- `title`:

  graph title

#### Returns

a ggplot object

------------------------------------------------------------------------

### `CrossValidation$clone()`

The objects of this class are cloneable with this method.

#### Usage

    CrossValidation$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
