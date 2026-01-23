/*
 * Author: Julien CHIQUET, INRAE
 *         julien.chiquet@inrae.fr
 *         Statistique et Génome
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_REGRESSION_DATA_H
#define _quadrupen_REGRESSION_DATA_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// TODO: use template to handle dense or sparse encoding (mat/sp_mat in armadillo)
class RegressionData {
  friend class GenericRegularizer ;
  friend class RidgeRegression ;
  friend class BoundedRegression ;
  friend class Optimizer ;
  // friend class ACTIVE_DATA; // NEED ACCES TO x, xt, y, xty, n, p

private:
  // DATA VARIABLES FOR REGRESSION PURPOSE
  uword  n_        ; // sample size
  uword  p_        ; // # of features
  mat    X_        ; // matrix of predictors
  vec    y_        ; // vector of response
  sp_mat S_        ; // Structuring matrix
  vec    weights_  ; // observation weights
  bool   centered_ ; // should intercept be considered?
  bool   scaled_   ; // should predictors be standardized?
  vec    XTy_      ; // responses to predictors vector
  vec    X_bar_    ; // mean of the predictors
  vec    norm_X_   ; // norm of the predictors
  double y_bar_    ; // mean of the response
  double norm_y_   ; // norm of the response
  
public:

  RegressionData(const Environment& dataModel, const bool& center, const bool& scale);
  
  // DATA NORMALIZATION
  // intercept treatment, predictor standardization, observation weights
  void standardize(const bool center, const bool scale);

  void scale_struct(const bool lambda);
  
  // getter functions to access private members
  const uword& n()          const { return n_        ; }
  const uword& p()          const { return p_        ; }
  const mat& X()            const { return X_        ; }
  const vec& XTy()          const { return XTy_      ; }
  const vec& y()            const { return y_        ; }
  const vec& norm_X()       const { return norm_X_   ; }
  const vec& X_bar()        const { return X_bar_    ; }
  const double& y_bar()     const { return y_bar_    ; }
  const double& norm_y()    const { return norm_y_   ; }
  const bool& is_centered() const { return centered_ ; }
  const bool& is_scaled()   const { return scaled_   ; }
  
};

#endif
