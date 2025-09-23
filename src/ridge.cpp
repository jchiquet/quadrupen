/*
 * Author: Julien CHIQUET
 *         Statistique et Génome*
 *         MIA PAris-Saclay
 */

#include "ridge.h"

using namespace Rcpp;
using namespace arma;

RIDGE::RIDGE(const arma::mat X, 
             const arma::vec Y,
             const arma::mat C, 
             const arma::vec WEIGHTS,
             const arma::uword VERBOSE) :
  x        (X), // predictors
  y        (Y), // responses
  c        (trimatu(C)), // Cholesky of L
  w        (WEIGHTS), // observation weights
  verbose  (VERBOSE)  // verbose mode (unused)
{}


// [[Rcpp::export]]
Rcpp::List ridge_cpp(
    const arma::mat& X  , // matrix of features
	  const arma::vec& Y        , // vector of response
	  const arma::mat& C        , // Cholesky decompostion of the structuring matrix
	  SEXP LAMBDA,
	  const arma::uword& NLAMBDA  ,
	  const double& LAMBDAMIN,
	  const double& LAMBDAMAX,
	  const bool& INTERCEPT,
	  const bool& NORMALIZE,
	  const arma::vec& WEIGHTS  ,
	  const arma::uword& VERBOSE  ) {

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
  mu = ybar - beta * xbar.t() ;
}

void RIDGE::standardize(const bool INTERCEPT, const bool NORMALIZE) {

  bool intercept  = (INTERCEPT) ; // boolean for intercept mode
  bool normalize  = (NORMALIZE) ; // boolean for standardizing the predictor
  uword n = x.n_rows;
  uword p = x.n_cols;

  if (intercept == 1) {
    xbar = mean(x, 0);
    ybar = mean(y) ;
    x.each_row() -= xbar ;
  } else {
    xbar = zeros<rowvec>(p) ;
    ybar = 0;
  }

  if (normalize == 1) {
    normx = sqrt(sum(square(x),0));
    x.each_row() /= normx ;
  } else {
    normx = ones<rowvec>(p);
  }
  normy = sqrt(sum(square(y))) ;

  if (intercept == 1) {
    xty = trans(y-ybar) * x ;
  } else {
    xty = y.t()*x ;
  }
}

void RIDGE::get_lambda(SEXP LAMBDA, 
                       const arma::uword NLAMBDA,
                       const double LAMBDAMIN,
                       const double LAMBDAMAX) {

  if (LAMBDA != R_NilValue) {
    lambda  = as<vec>(LAMBDA)  ;
  } else {
    lambda = exp10(linspace(log10(LAMBDAMIN), log10(LAMBDAMAX), NLAMBDA)) ;
  }

}
