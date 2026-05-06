# Data Class

Class for storing data and various fixed quantity

## Public fields

- `X`:

  matrix of regressor

- `y`:

  vector of response

- `C_inv`:

  Inverse of the Cholesky decomposition of S

- `S`:

  SDP structuring matrix

- `wy`:

  vector of observation weights

## Active bindings

- `d`:

  number of regressor

- `n`:

  sample size

- `sparse_encoding`:

  logical indicating if the matrix of regressor is sparsely encoded

- `varnames`:

  character, the names of the covariates/regressors

- `normx`:

  norm of each column of X

## Methods

### Public methods

- [`DataModel$new()`](#method-DataModel-initialize)

- [`DataModel$CholStruct()`](#method-DataModel-CholStruct)

- [`DataModel$splitTrainTest()`](#method-DataModel-splitTrainTest)

- [`DataModel$splitSubSamples()`](#method-DataModel-splitSubSamples)

- [`DataModel$clone()`](#method-DataModel-clone)

------------------------------------------------------------------------

### `DataModel$new()`

constructor for DataModel

#### Usage

    DataModel$new(
      covariates,
      outcome,
      cov_struct,
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

- `obs_weights`:

  vector of observations weights

- `check_args`:

  logical, should args be check at initialization?

- `cov_weights`:

  vector of covariates/regressors weights

------------------------------------------------------------------------

### `DataModel$CholStruct()`

Compute Cholesky factorization of the Structuring matrix

#### Usage

    DataModel$CholStruct()

------------------------------------------------------------------------

### `DataModel$splitTrainTest()`

a function splitting the data into train and test folds

#### Usage

    DataModel$splitTrainTest(
      nfolds = 10,
      folds = split(sample(1:self$n), rep(1:nfolds, length = self$n))
    )

#### Arguments

- `nfolds`:

  the number of folds

- `folds`:

  a list of vectors describing the folds (optional)

#### Returns

a list with train and test data and id.

------------------------------------------------------------------------

### `DataModel$splitSubSamples()`

a function splitting data into subsamples

#### Usage

    DataModel$splitSubSamples(
      n_subsamples = 50,
      subsample_size = floor(self$n/2),
      subsamples = replicate(n_subsamples, sample(1:self$n, subsample_size), simplify =
        FALSE),
      weakness = 1
    )

#### Arguments

- `n_subsamples`:

  the number of subsamples

- `subsample_size`:

  the subsample size

- `subsamples`:

  list with vector of subsamples (optional)

- `weakness`:

  coefficient for randomly weighting the regressor, default to 1

#### Returns

a list of DataModel, resampling of the original

------------------------------------------------------------------------

### `DataModel$clone()`

The objects of this class are cloneable with this method.

#### Usage

    DataModel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
