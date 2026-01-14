# Package index

## Main functions for user

Functions for fitting various structured regularization models

- [`bounded.reg()`](https://jchiquet.github.io/quadrupen/reference/bounded.reg.md)
  : Fit a linear model with infinity-norm plus ridge-like regularization
- [`elastic.net()`](https://jchiquet.github.io/quadrupen/reference/elastic.net.md)
  : Fit a linear model with elastic-net regularization
- [`fusedlasso()`](https://jchiquet.github.io/quadrupen/reference/fusedlasso.md)
  : A function for fitting generalized fused-Lasso problems
- [`lasso()`](https://jchiquet.github.io/quadrupen/reference/lasso.md) :
  Fit a linear model with lasso regularization
- [`ridge()`](https://jchiquet.github.io/quadrupen/reference/ridge.md) :
  Fit a linear model with a structured ridge regularization

## Classes and methods for handling the Quadrupen fits

R6 Classes for the user to manipulate the ouput of the main functions

- [`QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
  : Class "QuadrupenFit"
- [`ElasticNetFit`](https://jchiquet.github.io/quadrupen/reference/ElasticNetFit.md)
  : Class "ElasticNetFit"
- [`FusedLassoFit`](https://jchiquet.github.io/quadrupen/reference/FusedLassoFit.md)
  : Class "FusedLassoFit"
- [`BoundedRegressionFit`](https://jchiquet.github.io/quadrupen/reference/BoundedRegressionFit.md)
  : Class "BoundedRegression"
- [`RidgeRegressionFit`](https://jchiquet.github.io/quadrupen/reference/RidgeRegressionFit.md)
  : Class "RidgeRegressionFit"

## Auxiliaries Classes and methods

R6 Classes used as fields of the main Class Quadrupen fit

- [`DataModel`](https://jchiquet.github.io/quadrupen/reference/DataModel.md)
  : Data Class
- [`CrossValidation`](https://jchiquet.github.io/quadrupen/reference/CrossValidation.md)
  : Class CrossValidation
- [`StabilityPath`](https://jchiquet.github.io/quadrupen/reference/StabilityPath.md)
  : Class StabilityPath
- [`InformationCriteria`](https://jchiquet.github.io/quadrupen/reference/InformationCriteria.md)
  : Class InformationCriteria

## Auxiliaries S3 methods

S3 methods, aliases for the most useful R6 methods.

- [`criteria()`](https://jchiquet.github.io/quadrupen/reference/criteria.md)
  : Penalized criteria based on estimation of degrees of freedom
- [`cross_validate()`](https://jchiquet.github.io/quadrupen/reference/cross_validate.md)
  : Cross-validation for Quadrupen object
- [`stability()`](https://jchiquet.github.io/quadrupen/reference/stability.md)
  : Stability selection for Quadrupen object
- [`residuals(`*`<QuadrupenFit>`*`)`](https://jchiquet.github.io/quadrupen/reference/residuals.QuadrupenFit.md)
  : Extract model residuals
- [`coef(`*`<QuadrupenFit>`*`)`](https://jchiquet.github.io/quadrupen/reference/coef.QuadrupenFit.md)
  : Extract model coefficients
- [`predict(`*`<QuadrupenFit>`*`)`](https://jchiquet.github.io/quadrupen/reference/predict.QuadrupenFit.md)
  : Perform model prediction
- [`fitted(`*`<QuadrupenFit>`*`)`](https://jchiquet.github.io/quadrupen/reference/fitted.QuadrupenFit.md)
  : Extracts model fitted values
- [`deviance(`*`<QuadrupenFit>`*`)`](https://jchiquet.github.io/quadrupen/reference/deviance.QuadrupenFit.md)
  : Extract model deviance
- [`isQuadrupenFit()`](https://jchiquet.github.io/quadrupen/reference/isQuadrupenFit.md)
  : Auxiliary functions to check the given class of an object
