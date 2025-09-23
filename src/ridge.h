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

class RIDGE {
 public:
  RIDGE(const arma::mat X,
        const arma::vec Y,
        const arma::mat C,
        const arma::vec WEIGHTS,
        const uword VERBOSE);

  void standardize(const bool INTERCEPT, const bool NORMALIZE);

  void get_lambda(SEXP LAMBDA, 
                  const arma::uword NLAMBDA,
                  const double LAMBDAMIN,
                  const double LAMBDAMAX);

  void LetsRoll();

  // various function to access private members
  const arma::mat & get_coef() const { return beta; }
  const arma::vec & get_lambda() const { return lambda; }
  const arma::vec & get_df() const { return df; }
  const arma::rowvec & get_normx() const { return normx; }
  const arma::vec & get_mu() const { return mu; }

 private:
  arma::mat x       ; // matrix of predictors
  arma::vec y       ; // vector of reponses
  arma::mat c       ; // Cholesky decomposition of the structuring matrix
  arma::mat w       ; // vector of obervation weights
  int verbose       ; // integer for verbose mode

  arma::rowvec xty   ; // reponses to predictors vector
  arma::rowvec xbar  ; // mean of the predictors
  arma::rowvec normx ; // norm of the predictors
  double normy       ; // norm of the response
  double ybar        ; // mean of the response

  arma::vec lambda  ; // vector of ridge tuning parameters
  arma::mat beta    ; // matrix of coefficients
  arma::vec df      ; // degrees of freedom along the path
  arma::vec mu      ; // vector of intercept term

};

#endif

