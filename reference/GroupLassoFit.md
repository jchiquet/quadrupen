# Class "GroupLassoFit"

Class "GroupLassoFit"

Class "GroupLassoFit"

## Details

Class of object returned by the fitting function
[`group.lasso()`](https://jchiquet.github.io/quadrupen/reference/group.lasso.md).
Inherits fields and methods of
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)

## See also

[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
[`group.lasso()`](https://jchiquet.github.io/quadrupen/reference/group.lasso.md)

## Super class

[`quadrupen::QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
-\> `GroupLassoFit`

## Active bindings

- `penalty`:

  character describing the regularizer/penalty

- `group`:

  vector of integers indicating group belonging

- `type`:

  string indicating whether the \\\ell_1/\ell_2\\ or the
  \\\ell_1/\ell\_\infty\\ group-Lasso must be fitted.

## Methods

### Public methods

- [`GroupLassoFit$new()`](#method-GroupLassoFit-new)

- [`GroupLassoFit$clone()`](#method-GroupLassoFit-clone)

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

Initialize a
[`ElasticNetFit`](https://jchiquet.github.io/quadrupen/reference/ElasticNetFit.md)
model

#### Usage

    GroupLassoFit$new(data, intercept, group, type, regParam)

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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    GroupLassoFit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
