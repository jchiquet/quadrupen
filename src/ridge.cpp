/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "quadrupen_headers.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List ridge_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // boolean for intercept
    const List        &regParam    , // config of the optimisation 
    const List        &controlFit  // config of the optimisation 
) {

  const uword n             = dataModel["n"]  ; // sample size
  const uword p             = dataModel["d"]  ; // problem size
  const SEXP &R_X           = dataModel["X"]  ; // design matrix
  const arma::vec &y        = dataModel["y"]  ; // response vector
  const arma::mat& C        = dataModel["C"]      ; // Cholesky decompostion of the structuring matrix
  const arma::vec &wobs     = dataModel["wy"] ; // observation weights (not use at the moment)
  const bool sparse         = dataModel["sparse_encoding"] ; // boolean for sparse mode
  
  const SEXP R_LAMBDA       = regParam["l2"]         ; // vector of L1 penalties ; if NULL, automatically set
  const arma::vec &wlambda  = regParam["l2_weights"] ; // l1-penalty weights
  uword n_lambda            = regParam["n_lambda"]   ; // # of l1-penalty levels
  const double min_ratio    = regParam["min_ratio"]  ; // minimum penlaty value as a ratio of lambda1 max
  const double lambda_max   = regParam["lambda_max"]  ; // minimum penlaty value as a ratio of lambda1 max
  const bool normalize      = controlFit["normalize"] ; // boolean for standardizing the predictor

  vec    xty   ; // responses to predictors vector
  vec    xbar  ; // mean of the predictors
  vec    normx ; // norm of the predictors
  double normy ; // norm of the response
  double ybar  ; // mean of the response

  mat x        ;
  mat xt       ;
  x = as<mat>(R_X) ;
  standardize(x, y, intercept, normalize, wlambda, xty, normx, normy, xbar, ybar) ;

  // Vector of regularization
  vec lambda = get_lambda1(R_LAMBDA, n_lambda, min_ratio, lambda_max) ;
  n_lambda = lambda.n_elem ; // # of penalty levels
  
  // ==============================================================
  // INSTANTIATE THE RIDGE PATH
  arma::mat beta     ; // matrix of coefficients
  arma::vec df       ; // degrees of freedom along the path
  arma::vec mu       ; // vector of intercept term
  
  // ==============================================================
  // COMPUTE THE PATH OF SOLUTIONS
  arma::mat cinv = inv(trimatu(C)) ; // inverting the Cholesky decomp. of the structuring matrix
  
  // SVD DECOMPOSITION OF ( X * C^-1)
  arma::vec eta ; // eigen value of X cinv
  arma::mat U   ; // left singular vectors of X
  arma::mat V   ; // right singular vectors of X
  svd_econ(U, eta, V, (x.each_row() - xbar.t())*cinv) ;
  arma::mat cinvV = cinv * V ;
  arma::mat Uty = trans(U) * y ;
  
  beta.resize(lambda.n_elem, x.n_cols);
  df.resize(lambda.n_elem, 1);
  
  for (uword i=0; i<lambda.n_elem; i++) {
    // computing the structured ridge estimate
    beta.row(i) = trans( (cinvV * diagmat(eta/(square(eta) + lambda(i))) * Uty) / (normx % wlambda ));
    // computing the estimated degrees of freedom
    df(i) = sum(square(eta)/(square(eta) + lambda(i)));
  }
  
  // estimating the intercept term
  mu = ybar - beta * xbar ;

  return List::create(
    Named("tuning_param") = List::create(
      Named("l2") = lambda,
      Named("l1") = 0 
    ),
    Named("beta") = beta,
    Named("mu")   = mu  ,
    Named("df")   = df  ,
    Named("monitoring")  = List::create()
  );
  
}
