# Perform model prediction

Predict response for new sample based on the current model

## Usage

``` r
# S3 method for class 'QuadrupenFit'
predict(object, newx = NULL, selection = NULL, ...)
```

## Arguments

- object:

  an R object to evaluate

- newx:

  matrix of new values for the regressor with which to predict. If
  omitted, the fitted values are used.

- selection:

  either a character (model selection criteria) of a scalar (lambda
  value)

- ...:

  not used, only here for S3 compatibility

## Value

a vector of predicted value
