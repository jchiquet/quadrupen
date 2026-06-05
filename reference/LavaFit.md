# Class "LavaFit"

Class of object returned by the fitting function
[`lava()`](https://jchiquet.github.io/quadrupen/reference/lava.md).
Inherits fields and methods of
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

## See also

[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
[`lava()`](https://jchiquet.github.io/quadrupen/reference/lava.md)

## Super class

[`QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
-\> `LavaFit`

## Active bindings

- `penalty`:

  character describing the regularizer/penalty

- `lambda1`:

  vector of tuning parameters for the l1 penalty (sparse component)

- `lambda2`:

  vector of tuning parameters for the l2 penalty (dense component)

- `sparse_coef`:

  sparse part of the decomposition of the coefficients

- `dense_coef`:

  dense part of the decomposition of the coefficients

- `debias`:

  logical, should we rely on the debias coefficient of the regularizer
  (if available) or not

## Methods

### Public methods

- [`LavaFit$new()`](#method-LavaFit-initialize)

- [`LavaFit$fit()`](#method-LavaFit-fit)

- [`LavaFit$plot_path()`](#method-LavaFit-plot_path)

- [`LavaFit$clone()`](#method-LavaFit-clone)

Inherited methods

- [`QuadrupenFit$criteria()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-criteria)
- [`QuadrupenFit$cross_validate()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-cross_validate)
- [`QuadrupenFit$get_model()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-get_model)
- [`QuadrupenFit$plot()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-plot)
- [`QuadrupenFit$predict()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-predict)
- [`QuadrupenFit$print()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-print)
- [`QuadrupenFit$show()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-show)
- [`QuadrupenFit$stability()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-stability)

------------------------------------------------------------------------

### `LavaFit$new()`

Initialize a `LavaFit` model

#### Usage

    LavaFit$new(data, intercept, regParam)

#### Arguments

- `data`:

  a
  [`DataModel`](https://jchiquet.github.io/quadrupen/reference/DataModel.md)
  object

- `intercept`:

  a logical; should an intercept be included in the mode?

- `regParam`:

  a list with two elements, a vector and a scalar, for the
  regularization

------------------------------------------------------------------------

### `LavaFit$fit()`

function performing the optimization

#### Usage

    LavaFit$fit(control)

#### Arguments

- `control`:

  list controlling the optimization process Plot method for lava
  regularization path

------------------------------------------------------------------------

### `LavaFit$plot_path()`

Produce a plot of the solution path of a LavaFit object.

#### Usage

    LavaFit$plot_path(
      xvar = c("lambda", "fraction", "df"),
      log_scale = TRUE,
      component = "both",
      title = paste("Lava path:", component, "component(s)"),
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

- `component`:

  a character indicating the component to plot: both (sum of sparse and
  dense), sparse or dense. Default to both.

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

------------------------------------------------------------------------

### `LavaFit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    LavaFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
