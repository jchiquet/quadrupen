# Package index

## Main functions for user

Functions for fitting various structured regularization models

- [`bounded_reg()`](https://jchiquet.github.io/quadrupen/reference/bounded_reg.md)
  [`bounded.reg()`](https://jchiquet.github.io/quadrupen/reference/bounded_reg.md)
  : Fit a linear model with infinity-norm plus ridge-like regularization
- [`sparse_lm()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md)
  [`elastic.net()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md)
  [`elastic_net()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md)
  [`lasso()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md)
  [`mcp()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md)
  [`scad()`](https://jchiquet.github.io/quadrupen/reference/sparse_lm.md)
  : Fit a linear model with sparse regularization
- [`fused_lasso()`](https://jchiquet.github.io/quadrupen/reference/fused_lasso.md)
  : A function for fitting generalized fused-Lasso problems
- [`group_sparse_lm()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)
  [`group_lasso()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)
  [`group_l1linf()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)
  [`coop_lasso()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)
  [`sparse_group_lasso()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)
  [`sparse_group_l1linf()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)
  [`sparse_coop_lasso()`](https://jchiquet.github.io/quadrupen/reference/group_sparse_lm.md)
  : Fit a linear model with (sparse) group regularisation
- [`group_lava()`](https://jchiquet.github.io/quadrupen/reference/group_lava.md)
  : Fit a linear model with group-lava regularization
- [`lava()`](https://jchiquet.github.io/quadrupen/reference/lava.md) :
  Fit a linear model with lava regularization
- [`ridge()`](https://jchiquet.github.io/quadrupen/reference/ridge.md) :
  Fit a linear model with a structured ridge regularization

## Classes and methods for handling the Quadrupen fits

R6 Classes for the user to manipulate the ouput of the main functions

- [`QuadrupenFit`](https://jchiquet.github.io/quadrupen/reference/QuadrupenFit.md)
  : Class "QuadrupenFit"
- [`FusedLassoFit`](https://jchiquet.github.io/quadrupen/reference/FusedLassoFit.md)
  : Class "FusedLassoFit"
- [`SparseFit`](https://jchiquet.github.io/quadrupen/reference/SparseFit.md)
  : Class "SparseFit"
- [`SparseGroupFit`](https://jchiquet.github.io/quadrupen/reference/SparseGroupFit.md)
  : Class "SparseGroupFit"
- [`GroupLavaFit`](https://jchiquet.github.io/quadrupen/reference/GroupLavaFit.md)
  : Class "GroupLavaFit"
- [`LavaFit`](https://jchiquet.github.io/quadrupen/reference/LavaFit.md)
  : Class "LavaFit"
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
