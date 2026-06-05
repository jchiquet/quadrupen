# Class "GroupLavaFit"

Class of object returned by the fitting function
[`group_lava()`](https://jchiquet.github.io/quadrupen/reference/group_lava.md).
Inherits fields and methods of
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
and [LavaFit](https://jchiquet.github.io/quadrupen/reference/LavaFit.md)

## See also

[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
[`group_lava()`](https://jchiquet.github.io/quadrupen/reference/group_lava.md)

## Super classes

[`QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
-\>
[`LavaFit`](https://jchiquet.github.io/quadrupen/reference/LavaFit.md)
-\> `GroupLavaFit`

## Active bindings

- `lambda1`:

  vector of tuning parameters for the l1 group penalty (sparse
  component)

- `lambda2`:

  vector of tuning parameters for the l2 penalty (dense component)

- `penalty`:

  character describing the regularizer/penalty

- `group`:

  vector of integers indicating group belonging

- `type`:

  string indicating whether the \\\ell_1/\ell_2\\ or the
  \\\ell_1/\ell\_\infty\\ group-Lasso must be fitted.

## Methods

### Public methods

- [`GroupLavaFit$new()`](#method-GroupLavaFit-initialize)

- [`GroupLavaFit$clone()`](#method-GroupLavaFit-clone)

Inherited methods

- [`QuadrupenFit$criteria()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-criteria)
- [`QuadrupenFit$cross_validate()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-cross_validate)
- [`QuadrupenFit$get_model()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-get_model)
- [`QuadrupenFit$plot()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-plot)
- [`QuadrupenFit$predict()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-predict)
- [`QuadrupenFit$print()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-print)
- [`QuadrupenFit$show()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-show)
- [`QuadrupenFit$stability()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-stability)
- [`LavaFit$fit()`](https://jchiquet.github.io/quadrupen/reference/LavaFit.html#method-fit)
- [`LavaFit$plot_path()`](https://jchiquet.github.io/quadrupen/reference/LavaFit.html#method-plot_path)

------------------------------------------------------------------------

### `GroupLavaFit$new()`

Initialize a `GroupLavaFit` model

#### Usage

    GroupLavaFit$new(data, intercept, group, type, regParam)

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

### `GroupLavaFit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    GroupLavaFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
