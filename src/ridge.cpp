/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#include "ridge.h"

using namespace Rcpp;
using namespace arma;

RIDGE::RIDGE(SEXP X, SEXP Y, SEXP C, SEXP WEIGHTS, SEXP VERBOSE) :
  x        (as<mat>         (X)     ) , // predictors
  y        (as<vec>         (Y)     ) , // responses
  c        (trimatu(as<mat> (C)     )), // cholesky of L
  w        (as<vec>         (WEIGHTS)), // observation weights
  verbose  (as<int>         (VERBOSE))  // verbose mode (unused)
{}


SEXP ridge_cpp(SEXP X        , // matrix of features
	       SEXP Y        , // vector of response
	       SEXP C        , // Cholesky decompostion of the structuring matrix
	       SEXP LAMBDA   ,
	       SEXP NLAMBDA  ,
	       SEXP LAMBDAMIN,
	       SEXP LAMBDAMAX,
	       SEXP INTERCEPT,
	       SEXP NORMALIZE,
	       SEXP WEIGHTS  ,
	       SEXP VERBOSE  ) {

  // ==============================================================
  // INSTANTIATE THE RIDGE PATH
  RIDGE ridge(X, Y, C, WEIGHTS, VERBOSE);

  // ==============================================================
  // DATA NORMALIZATION
  ridge.standardize(INTERCEPT, NORMALIZE) ;

  // ==============================================================
  // GET LAMBDA VALUES
  ridge.get_lambda(LAMBDA, NLAMBDA, LAMBDAMIN, LAMBDAMAX) ;

  // ==============================================================
  // COMPUTE THE PATH OF SOLUTIONS
  ridge.LetsRoll();

  return List::create(Named("coefficients") = ridge.get_coef()  ,
		      Named("mu")           = ridge.get_mu()    ,
		      Named("normx")        = ridge.get_normx() ,
		      Named("lambda2")      = ridge.get_lambda(),
		      Named("df")           = ridge.get_df()    );

}

void RIDGE::LetsRoll() {

  arma::mat cinv = inv(trimatu(c)) ; // inverting the Cholesky decomp. of the structuring matrix

  // SVD DECOMPOSITION OF ( X * C^-1)
  arma::vec eta ; // eigen value of X cinv
  arma::mat U   ; // left singular vectors of X
  arma::mat V   ; // right singular vectors of X
  svd_econ(U, eta, V, x*cinv) ;
  
  arma::mat cinvV = cinv * V ;
  arma::mat Uty = trans(U) * y ;

  beta.resize(lambda.n_elem, x.n_cols);
  df.resize(lambda.n_elem, 1);

  for (int i; i<lambda.n_elem; i++) {
    // computing the structured ridge estimate
    beta.row(i) = trans(cinvV * diagmat(eta/(square(eta) + lambda(i))) * Uty) / normx;
    // computing the estimated degrees of freedom
    df(i) = sum(square(eta)/(square(eta) + lambda(i)));
  }

  // estimating the intercept term
  mu = ybar - beta * xbar;
}

void RIDGE::standardize(SEXP INTERCEPT, SEXP NORMALIZE) {

  bool intercept  = as<bool>   (INTERCEPT) ; // boolean for intercept mode
  bool normalize  = as<bool>   (NORMALIZE) ; // boolean for standardizing the predictor
  uword n = x.n_rows;
  uword p = x.n_cols;

  if (intercept == 1) {
    xbar = trans(rowvec(mean(x, 0)));
    ybar = mean(y) ;
    for (int i=0; i<p; i++) {
      x.col(i) -= xbar(i);
    }
  } else {
    xbar = zeros(p) ;
    ybar = 0;
  }

  if (normalize == 1) {
    normx = sqrt(trans(sum(square(x),0)));
    for (int i=0; i<p; i++) {
      x.col(i) /= normx(i);
    }
  } else {
    normx = ones(p);
  }
  normy = sqrt(sum(square(y))) ;

  if (intercept == 1) {
    xty = trans(trans(y-ybar) * x) ;
  } else {
    xty = trans(y.t()*x) ;
  }
}

void RIDGE::get_lambda(SEXP LAMBDA, SEXP NLAMBDA, SEXP LAMBDAMIN, SEXP LAMBDAMAX) {

  if (LAMBDA != R_NilValue) {
    lambda  = as<vec>(LAMBDA)  ;
  } else {
    lambda = exp10(linspace(log10(as<double>(LAMBDAMIN)), log10(as<double>(LAMBDAMAX)), as<uword>(NLAMBDA))) ;
  }

}
