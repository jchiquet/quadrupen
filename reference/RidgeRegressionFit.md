# Class "RidgeRegressionFit"

Class "RidgeRegressionFit"

Class "RidgeRegressionFit"

## Details

Class of object returned by the fitting function
[`ridge()`](https://jchiquet.github.io/quadrupen/reference/ridge.md).
Inherits fields and methods of
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md).

## See also

[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
[`ridge()`](https://jchiquet.github.io/quadrupen/reference/ridge.md)

## Super class

[`quadrupen::QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
-\> `RidgeRegressionFit`

## Active bindings

- `penalty`:

  character describing the regularizer/penalty

## Methods

### Public methods

- [`RidgeRegressionFit$new()`](#method-RidgeRegressionFit-new)

- [`RidgeRegressionFit$clone()`](#method-RidgeRegressionFit-clone)

Inherited methods

- [`quadrupen::QuadrupenFit$criteria()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-criteria)
- [`quadrupen::QuadrupenFit$cross_validate()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-cross_validate)
- [`quadrupen::QuadrupenFit$fit()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-fit)
- [`quadrupen::QuadrupenFit$get_model()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-get_model)
- [`quadrupen::QuadrupenFit$plot()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-plot)
- [`quadrupen::QuadrupenFit$plot_path()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-plot_path)
- [`quadrupen::QuadrupenFit$predict()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-predict)
- [`quadrupen::QuadrupenFit$print()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-print)
- [`quadrupen::QuadrupenFit$show()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-show)
- [`quadrupen::QuadrupenFit$stability()`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.html#method-stability)

------------------------------------------------------------------------

### Method `new()`

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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    RidgeRegressionFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
