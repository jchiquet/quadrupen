# Plot method for quadrupen objects

S3 plot methods for
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
[CrossValidation](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md)
and
[StabilityPath](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
objects, delegating to their respective R6 `$plot()` method.

## Usage

``` r
# S3 method for class 'QuadrupenFit'
plot(x, ...)

# S3 method for class 'CrossValidation'
plot(x, ...)

# S3 method for class 'StabilityPath'
plot(x, ...)
```

## Arguments

- x:

  a
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md),
  [CrossValidation](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md)
  or
  [StabilityPath](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
  object.

- ...:

  additional arguments passed to the underlying R6 `$plot()` method. For
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md):
  `type` (`"path"`, `"criteria"`, `"crossval"`, `"stability"`),
  `log_scale`, `labels`. For
  [CrossValidation](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md):
  `log_scale`, `title`. For
  [StabilityPath](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md):
  `xvar`, `title`, `labels`, `sel_mode`, `cutoff`, `PFER`, `nvarsel`.

## Value

a ggplot2 object.

## Functions

- `plot(CrossValidation)`: Plot method for a
  [CrossValidation](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md)
  object

- `plot(StabilityPath)`: Plot method for a
  [StabilityPath](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
  object
