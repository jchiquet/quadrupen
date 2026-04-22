# Changelog

## quadrupen 0.3-0 (2013-11)

- major updates
  - added structured ridge regression (regression penalized by
    structured l2 norm)
  - added corresponding functionality for cross-validation and stability
    path
- Minor:
  - removed the not very useful ‘reverse’ argument for plotting  
  - added demo for lasso, ridge and use of the criteria function

## quadrupen 0.2-4 (2014-01-16)

CRAN release: 2014-01-16

``` R
Minor:
+   memory leak corrected (sp_mat declaration)
+   linking to Rcpp/RcppArmadillo headers (requires R 3.0-2)
```

## quadrupen 0.2-4 (2013-11-01)

CRAN release: 2014-01-16

- added a ‘lasso’ function, simple wrapper to the elastic-net ‘function’
- added computation of degrees of freedom (for elastic net and bounded
  regression)
- added a method to compute penalized criteria (BIC/AIC) of a quadrupen
  fit, with plot

## quadrupen 0.2-2 (2013-04-08)

CRAN release: 2013-04-08

- minor fix to comply with recent ggplot2 updates.

## quadrupen 0.2-1 (2013-02-27)

CRAN release: 2013-02-27

- minor fix to pass CRAN check on Windows operating systems.

quadrupen 0.2-0 (2013-02-26)

- Major updates
  - added bounded regression (regression penalized by infinity norm +
    structured l2 norm)
  - added corresponding features for cross-validation and stability path
- Minor updates
- corrected wrong annotations of the stability path (PFER)
- handled normalization internally (‘normalize’ is no longer a
  parameter)
- more simple internal handling of penscales and correction of the
  rescaling of the intercept
- better use of multicore features
- handled runtime error exception in RcppArmadillo when the system is
  singular (end of the solution path) A consequence is quadrupen is less
  likely to crash due to user’s “bad” parametrisation
- simplification of the C++ code, bugs corrected, probably new ones
  added :-’(
- added ‘examples’ and ‘tests’ directories

## quadrupen 0.1-0 (2012-10-09)

CRAN release: 2012-10-10

- first build: structured elastic-net with (weighted) quadratic loss,
  cross-validation and stability selection methods.
