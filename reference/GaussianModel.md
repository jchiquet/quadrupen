# Gaussian Data Class

Class for storing data and various fixed quantity

## Super class

[`quadrupen::DataModel`](https://jchiquet.github.io/quadrupen/reference/DataModel.md)
-\> `GaussianModel`

## Public fields

- `xty`:

  sufficient statistics for the Gaussian Data

## Active bindings

- `name`:

  Type of data

- `rss`:

  residuals sum of squares of the responses sum((y - mean(y))^2)

## Methods

### Public methods

- [`GaussianModel$getSufficientStat()`](#method-GaussianModel-getSufficientStat)

- [`GaussianModel$show()`](#method-GaussianModel-show)

- [`GaussianModel$print()`](#method-GaussianModel-print)

- [`GaussianModel$splitTrainTest()`](#method-GaussianModel-splitTrainTest)

- [`GaussianModel$splitSubSamples()`](#method-GaussianModel-splitSubSamples)

- [`GaussianModel$clone()`](#method-GaussianModel-clone)

Inherited methods

- [`quadrupen::DataModel$CholStruct()`](https://jchiquet.github.io/quadrupen/reference/DataModel.html#method-CholStruct)
- [`quadrupen::DataModel$initialize()`](https://jchiquet.github.io/quadrupen/reference/DataModel.html#method-initialize)
- [`quadrupen::DataModel$standardize()`](https://jchiquet.github.io/quadrupen/reference/DataModel.html#method-standardize)

------------------------------------------------------------------------

### Method `getSufficientStat()`

compute sufficient statistics for the Gaussian Data

#### Usage

    GaussianModel$getSufficientStat()

------------------------------------------------------------------------

### Method `show()`

a print method

#### Usage

    GaussianModel$show()

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

User friendly print method

#### Usage

    GaussianModel$print()

------------------------------------------------------------------------

### Method `splitTrainTest()`

a function splitting the data into train and test folds

#### Usage

    GaussianModel$splitTrainTest(
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

### Method `splitSubSamples()`

a function splitting data into subsamples

#### Usage

    GaussianModel$splitSubSamples(
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

  coefficient for randmonly reweighting the regressor, default to 1

#### Returns

a list of DataModel, resampling of the original

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    GaussianModel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
