# Class "SparseGroupFit"

Class of object returned by the fitting function
[`group_sparse_lm()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md).
Inherits fields and methods of
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

## See also

[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
[`group_sparse_lm()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)

## Super class

[`QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
-\> `SparseGroupFit`

## Active bindings

- `penalty`:

  character describing the regularizer/penalty

- `group`:

  vector of integers indicating group belonging

- `type`:

  string the type of group-wise regularization applied

- `mixture_tuning`:

  mixture coefficient of the sparse group penalty

- `is_group_sparse`:

  boolean indicating if sparse group or group penalty is applied

## Methods

### Public methods

- [`SparseGroupFit$new()`](#method-SparseGroupFit-initialize)

- [`SparseGroupFit$clone()`](#method-SparseGroupFit-clone)

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

### `SparseGroupFit$new()`

Initialize a `SparseGroupFit` model

#### Usage

    SparseGroupFit$new(data, intercept, group, type, regParam)

#### Arguments

- `data`:

  a
  [`DataModel`](https://jchiquet.github.io/quadrupen/reference/DataModel.md)
  object

- `intercept`:

  a logical; should an intercept be included in the mode?

- `group`:

  vector of integers indicating group belonging.

- `type`:

  string indicating whether the \\\ell_1/\ell_2\\ or the
  \\\ell_1/\ell\_\infty\\ group-Lasso must be fitted.

- `regParam`:

  a list with two elements, a vector and a scalar, for the
  regularization

------------------------------------------------------------------------

### `SparseGroupFit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SparseGroupFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
