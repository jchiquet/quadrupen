#include "data_reg.hpp"

using namespace Rcpp;
using namespace arma;

REGRESSION_DATA::
REGRESSION_DATA(SEXP X          , // matrix of features
		SEXP Y          , // vector of response
		SEXP INTERCEPT  , // should intercept be considered?
		SEXP NORMALIZE  , // should predictors be standardized?
		SEXP PREWEIGHTS , // prediction weights
		SEXP OBSWEIGHTS): // maximal number of features allowed to enter the path

  // BASICALLY READING THE INPUT ARGUMENTS FROM R
  x          (as<mat>    (X))         ,
  y          (as<vec>    (Y))         ,
  intercept  (as<bool>   (INTERCEPT)) ,
  normalize  (as<bool>   (NORMALIZE)) ,
  preweights (as<vec>    (PREWEIGHTS)),
  obsweights (as<vec>    (OBSWEIGHTS)) {
  n  = x.n_rows ;
  p  = x.n_cols ;
}

void REGRESSION_DATA::standardize() {

  // TODO: OBSERVATION WEIGHTS
  
  if (intercept == 1) {
    xbar = mean(x).t();
    ybar = mean(y) ;
    for (int i=0; i<p; i++) {
      x.col(i) -= xbar(i);
    }
  } else {
    xbar = zeros(p) ;
    ybar = 0;
  }

  if (normalize == 1) {
    normx = sqrt(trans(sum(square(x))));
    for (int i=0; i<p; i++) {
      x.col(i) /= normx(i);
    }
  } else {
    normx = ones(p);
  }
  normy = sqrt(sum(square(y))) ;

  if (any(preweights != 1)) {
    x.each_row() /= preweights.t();
    xbar /= preweights;
  }

  xt = x.t();
  if (intercept == 1) {
    xty = xt * (y-ybar) ;
  } else {
    xty = xt * y ;
  }

}
