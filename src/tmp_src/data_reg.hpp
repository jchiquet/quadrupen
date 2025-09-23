/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_DATA_REG_H
#define _quadrupen_DATA_REG_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// class OPTIM_DATA ;

// TODO: use template to handle dense or sparse encoding (mat/sp_mat in armadillo)
class REGRESSION_DATA {
  friend class ACTIVE_DATA; // NEED ACCES TO x, xt, y, xty, n, p

private:
  // DATA VARIABLES FOR REGRESSION PURPOSE
  uword  n          ; // sample size
  uword  p          ; // # of features
  vec    xbar       ; // mean of the predictors
  vec    normx      ; // norm of the predictors
  double ybar       ; // mean of the response
  double normy      ; // norm of the response
  bool   intercept  ; // should intercept be considered?
  bool   normalize  ; // should predictors be standardized?
  vec    obsweights ; // observation weights
  vec    preweights ; // predictor weights
  mat    x          ; // matrix of predictors
  mat    xt         ; // transpose matrix of predictor once and keep it to save time
  vec    xty        ; // reponses to predictors vector
  vec    y          ; // vector of response

public:
  REGRESSION_DATA(SEXP X          , // matrix of features
		  SEXP Y          , // vector of response
		  SEXP INTERCEPT  , // should intercept be considered?
		  SEXP NORMALIZE  , // should predictors be standardized?
		  SEXP PREWEIGHTS , // prediction weights
		  SEXP OBSWEIGHTS); // observation weights

  // DATA NORMALIZATION
  // intercept treatment, predictor standardization, predictor weihgts and observation weights
  void standardize();

  // various function to acces private members
  const mat& X()   const { return x   ; }
  const vec& Y()   const { return y   ; }
  const vec& XTY() const { return xty ; }

};

#endif
