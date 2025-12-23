# Extract model residuals

Extracts model residuals from a
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
object

## Usage

``` r
# S3 method for class 'QuadrupenFit'
residuals(object, newx = NULL, newy, ...)
```

## Arguments

- object:

  a
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
  object

- newx:

  matrix of new values for the regressor with which to predict. If
  omitted, the fitted values are used.

- newy:

  vector of new values for the response with which to compute the
  residuals. If omitted, the fitted values are used.

- ...:

  not used, only here for S3 compatibility

## Value

Matrix of residuals, each column corresponding to a value of `lambda1`.
