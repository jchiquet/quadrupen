/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_RIDGE_H
#define _quadrupen_RIDGE_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "quadrupen_headers.hpp"

using namespace Rcpp;
using namespace arma;

RcppExport SEXP ridge_cpp(SEXP X        , // matrix of features
			  SEXP Y        , // vector of response
			  SEXP C        , // Cholesky decompostion of the structuring matrix
			  SEXP LAMBDA   ,
			  SEXP NLAMBDA  ,
			  SEXP LAMBDAMIN,
			  SEXP LAMBDAMAX,
			  SEXP INTERCEPT,
			  SEXP NORMALIZE,
			  SEXP WEIGHTS  ,
			  SEXP VERBOSE  ) ;

class RIDGE {
 public:
  RIDGE(SEXP X, SEXP Y, SEXP C, SEXP WEIGHTS, SEXP VERBOSE);

  void standardize(SEXP INTERCEPT, SEXP NORMALIZE);

  void get_lambda(SEXP LAMBDA, SEXP NLAMBDA, SEXP LAMBDAMIN, SEXP LAMBDAMAX);

  void LetsRoll();

  // various function to acces private members
  const arma::mat & get_coef() const { return beta; }
  const arma::vec & get_lambda() const { return lambda; }
  const arma::vec & get_df() const { return df; }
  const arma::vec & get_normx() const { return normx; }
  const arma::vec & get_mu() const { return mu; }

 private:
  arma::mat x       ; // matrix of predictors
  arma::vec y       ; // vector of reponses
  arma::mat c       ; // Cholesky decomposition of the structuring matrix
  arma::mat w       ; // vector of obervation weights
  int verbose       ; // integer for verbose mode

  arma::vec xty     ; // reponses to predictors vector
  arma::vec xbar    ; // mean of the predictors
  arma::vec normx   ; // norm of the predictors
  double normy      ; // norm of the response
  double ybar       ; // mean of the response

  arma::vec lambda  ; // vector of ridge tuning parameters
  arma::mat beta    ; // matrix of coefficients
  arma::vec df      ; // degrees of freedom along the path
  arma::vec mu      ; // vector of intercept term

};

#endif

