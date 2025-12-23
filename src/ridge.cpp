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
    const List        &tuningParam , // List of tuning parameters
    const List        &control      // config of the optimisation 
) {
  
  arma::vec lambda          = tuningParam["l2"]   ; // vector of L2 penalties  
  const arma::mat &X        = dataModel["X"]      ; // design matrix
  const arma::vec &y        = dataModel["y"]      ; // response vector
  const arma::mat& C        = dataModel["C"]      ; // Cholesky decompostion of the structuring matrix
  const arma::rowvec &xbar  = dataModel["mean_X"] ; // mean of the predictors
  const arma::rowvec &normx = dataModel["norm_X"] ; // norm of the predictors
  const double ybar         = dataModel["mean_y"] ; // mean of the predictors
  
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
  svd_econ(U, eta, V, (X.each_row() - xbar)*cinv) ;
  arma::mat cinvV = cinv * V ;
  arma::mat Uty = trans(U) * y ;
  
  beta.resize(lambda.n_elem, X.n_cols);
  df.resize(lambda.n_elem, 1);
  
  for (uword i; i<lambda.n_elem; i++) {
    // computing the structured ridge estimate
    beta.row(i) = trans(cinvV * diagmat(eta/(square(eta) + lambda(i))) * Uty) / normx;
    // computing the estimated degrees of freedom
    df(i) = sum(square(eta)/(square(eta) + lambda(i)));
  }
  
  // estimating the intercept term
  mu = ybar - beta * xbar.t() ;

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
