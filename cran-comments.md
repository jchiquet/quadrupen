## quadrupen 1.1-0	(2026-09-17)

- performance: faster active set algorithms (no reallocation of the Gram matrices, several
  variables activated at once, restarted FISTA, exact block solver for group models, no p x p
  algebra for ridge and lava with a diagonal structure); typical speed-ups over 1.0-0 range
  from x2 to x36 depending on the model and the size of the active set
- bug fixes: bounded regression with the default solver could return suboptimal solutions;
  wrong proximal operator for the l1/linf group penalty; robustness when the active set
  exceeds the rank of the design
- new control option `maxadd`

## Tested environments

* tested locally on Ubuntu Linux 24.04 LTS, R 4.6.1, GCC 13.3

* tested remotely with github-action, all status OK

- Linux ubuntu 24.04, R-release
- Linux ubuntu 24.04, R-oldrel
- Linux ubuntu 24.04, R-devel
- Windows Server 2025, R-release, 64 bit
- macOS 15, R-release

* tested remotely with win-builder (R-release, R-devel, R-oldrelease), all status OK

* additionally tested remotely with R-hub v2 (memory-checking platforms and alternate
  BLAS/LAPACK, given the amount of C++/Rcpp/RcppArmadillo code)

- clang-ASAN: OK
- clang-UBSAN: OK
- gcc-ASAN: OK
- ATLAS (alternate BLAS/LAPACK implementation): OK
- valgrind: the example suite completed cleanly (0 errors, 0 bytes leaked, definitely/
  indirectly/possibly); the testthat suite is too slow under valgrind's instrumentation to
  finish within GitHub Actions' 6h job limit and was cancelled mid-run, with no error reported
  up to that point
- rchk: all findings are located inside Rcpp/RcppArmadillo header code (PROTECT/UNPROTECT
  static analysis of e.g. Armor.h, Shield.h, S4.h), none in the package's own sources; a
  known false-positive-prone pattern of this tool on Rcpp-based packages

## Local R CMD check results

── R CMD check results ── quadrupen 1.1-0 ────
(`R CMD check --as-cran --no-manual`, OMP_NUM_THREADS = OPENBLAS_NUM_THREADS = 2)

0 errors | 0 warnings | 1 note

* checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s): '-mno-omit-leaf-frame-pointer'

  This flag comes from the local R configuration (Ubuntu's Makeconf), not from the package.
