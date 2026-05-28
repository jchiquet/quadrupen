# Extract model residuals

Extracts model residuals from a
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
object

## Usage

``` r
# S3 method for class 'QuadrupenFit'
residuals(object, newx = NULL, newy = NULL, ...)
```

## Arguments

- object:

  a
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
  object

- newx:

  matrix of new covariates for out-of-sample residuals. Must be provided
  together with `newy`. If `NULL` (default), training residuals are
  returned.

- newy:

  vector of new responses for out-of-sample residuals. Must be provided
  together with `newx`. If `NULL` (default), training residuals are
  returned.

- ...:

  not used, only here for S3 compatibility

## Value

Matrix of residuals, each column corresponding to a value of `lambda1`.
