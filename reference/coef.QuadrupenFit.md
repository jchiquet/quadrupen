# Extract model coefficients

Extracts model coefficients from a
[QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
object

## Usage

``` r
# S3 method for class 'QuadrupenFit'
coef(object, selection = NULL, ...)
```

## Arguments

- object:

  a
  [QuadrupenFit](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
  object

- selection:

  either a character (model selection criteria) of a scalar (lambda
  value)

- ...:

  not used, only here for S3 compatibility

## Value

a vector of coefficients
