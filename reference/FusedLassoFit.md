# Class "FusedLassoFit"

Class of object returned by the fitting function
[`fusedlasso()`](https://jchiquet.github.io/quadrupen/reference/fusedlasso.md).
Inherits fields and methods of
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md).

## See also

[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
[`bounded.reg()`](https://jchiquet.github.io/quadrupen/reference/bounded.reg.md)

## Super class

[`QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
-\> `FusedLassoFit`

## Active bindings

- `penalty`:

  character describing the regularizer/penalty

## Methods

### Public methods

- [`FusedLassoFit$new()`](#method-FusedLassoFit-initialize)

- [`FusedLassoFit$clone()`](#method-FusedLassoFit-clone)

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

### `FusedLassoFit$new()`

Initialize a `FusedLassoFit` model

#### Usage

    FusedLassoFit$new(data, intercept, regParam)

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

### `FusedLassoFit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    FusedLassoFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
