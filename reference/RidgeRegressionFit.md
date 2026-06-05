# Class "RidgeRegressionFit"

Class of object returned by the fitting function
[`ridge()`](https://jchiquet.github.io/quadrupen/reference/ridge.md).
Inherits fields and methods of
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md).

## See also

[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
[`ridge()`](https://jchiquet.github.io/quadrupen/reference/ridge.md)

## Super class

[`QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
-\> `RidgeRegressionFit`

## Active bindings

- `penalty`:

  character describing the regularizer/penalty

- `lambda2`:

  vector of tuning parameters for the l2 penalty

## Methods

### Public methods

- [`RidgeRegressionFit$new()`](#method-RidgeRegressionFit-initialize)

- [`RidgeRegressionFit$clone()`](#method-RidgeRegressionFit-clone)

Inherited methods

- [`QuadrupenFit$criteria()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-criteria)
- [`QuadrupenFit$cross_validate()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-cross_validate)
- [`QuadrupenFit$fit()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-fit)
- [`QuadrupenFit$get_model()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-get_model)
- [`QuadrupenFit$plot()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-plot)
- [`QuadrupenFit$plot_path()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-plot_path)
- [`QuadrupenFit$predict()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-predict)
- [`QuadrupenFit$print()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-print)
- [`QuadrupenFit$show()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-show)
- [`QuadrupenFit$stability()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-stability)

------------------------------------------------------------------------

### `RidgeRegressionFit$new()`

Initialize a `RidgeRegressionFit` model

#### Usage

    RidgeRegressionFit$new(data, intercept, regParam)

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

### `RidgeRegressionFit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    RidgeRegressionFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
