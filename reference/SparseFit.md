# Class "SparseFit"

Class of object returned by the fitting function
[`elastic_net()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md).
Inherits fields and methods of
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

## See also

[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
[`elastic_net()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md)

## Super class

[`QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
-\> `SparseFit`

## Active bindings

- `lambda1`:

  vector of tuning parameters for the l1 penalty

- `lambda2`:

  vector of tuning parameters for the l2 penalty

- `penalty`:

  character describing the regularizer/penalty

- `type`:

  string the type of group-wise regularization applied

- `unbiasing_tuning`:

  unbiasing coefficient of the MCP or SCAD penalties

## Methods

### Public methods

- [`SparseFit$new()`](#method-SparseFit-initialize)

- [`SparseFit$clone()`](#method-SparseFit-clone)

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

### `SparseFit$new()`

Initialize a `SparseFit` model

#### Usage

    SparseFit$new(data, intercept, type, regParam)

#### Arguments

- `data`:

  a
  [`DataModel`](https://jchiquet.github.io/quadrupen/reference/DataModel.md)
  object

- `intercept`:

  a logical; should an intercept be included in the mode?

- `type`:

  string the type of group-wise regularization applied

- `regParam`:

  a list with two elements, a vector and a scalar, for the
  regularization

------------------------------------------------------------------------

### `SparseFit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SparseFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
