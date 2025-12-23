# Data Class

Class for storing data and various fixed quantity

## Public fields

- `X`:

  matrix of regressor

- `y`:

  vector of response

- `C`:

  Cholesky decomposition of the structring matrix S

- `S`:

  SDP structuring matrix

- `wx`:

  vector of regressor weights

- `wy`:

  vector of observation weights

- `mean_X`:

  vector of (normalized) regressor means

- `norm_X`:

  vector of regressor norms

- `mean_y`:

  mean of the response vector

- `norm_y`:

  norm of the response vector

## Active bindings

- `d`:

  number of regressor

- `n`:

  sample size

- `is_centered`:

  logical indicating if the data has been centered

- `is_scaled`:

  logical indicating if the data has been scaled

- `sparse_encoding`:

  logical indicating if the matrix of regressor is sparsely encoded

- `varnames`:

  character, the names of the covariates/regressors

## Methods

### Public methods

- [`DataModel$new()`](#method-DataModel-new)

- [`DataModel$standardize()`](#method-DataModel-standardize)

- [`DataModel$CholStruct()`](#method-DataModel-CholStruct)

- [`DataModel$clone()`](#method-DataModel-clone)

------------------------------------------------------------------------

### Method `new()`

constructor for DataModel

#### Usage

    DataModel$new(
      covariates,
      outcome,
      cov_struct,
      cov_weights = rep(1, ncol(covariates)),
      obs_weights = rep(1, length(outcome)),
      check_args = TRUE
    )

#### Arguments

- `covariates`:

  matrix of covariates/regressors

- `outcome`:

  vector of outcome/response

- `cov_struct`:

  sdp matrix structuring the covariates/regressors

- `cov_weights`:

  vector of covariates/regressors weights

- `obs_weights`:

  vector of observations weights

- `check_args`:

  logical, should args be check at initialization?

------------------------------------------------------------------------

### Method `standardize()`

Perform standardization of the data an store the auxiliaries quantities

#### Usage

    DataModel$standardize(intercept, normalize)

#### Arguments

- `intercept`:

  logical, is there an intercept in the model?

- `normalize`:

  logical, shall the regressor be standardized?

------------------------------------------------------------------------

### Method `CholStruct()`

Compute Cholesky factorization of the Structuring matrix

#### Usage

    DataModel$CholStruct()

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    DataModel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
