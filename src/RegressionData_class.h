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

template <typename matrix>  class GenericRegularizer ;
template <typename matrix>  class Optimizer ;

// Use template to handle dense or sparse encoding (mat/sp_mat in armadillo)
template <typename matrix>
class RegressionData {
  friend class GenericRegularizer<matrix> ;
  friend class RidgeRegression ;
  friend class BoundedRegression ;
  friend class Optimizer<matrix> ;

private:
  // DATA VARIABLES FOR REGRESSION PURPOSE
  uword  n_        ; // sample size
  uword  p_        ; // # of features
  matrix X_        ; // matrix of predictors
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

  void scale_struct(vec weights) {
    sp_mat diag_S = spdiags(weights, ivec({0}), p_, p_) ;
    S_ = diag_S * S_ * diag_S  ; 
  };
  
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

template <typename matrix>
RegressionData<matrix>::RegressionData(
  const Environment& dataModel,
  const bool& center, const bool& scale):
  n_       (as<uword>  (dataModel["n"])) , // sample size
  p_       (as<uword>  (dataModel["d"])) , // # of features
  y_       (as<vec>    (dataModel["y"])) , // response vector
  S_       (as<sp_mat> (dataModel["S"])) , // structuring matrix
  weights_ (as<vec>    (dataModel["wy"]))  // observation weights
{
  if (Rf_isS4(dataModel["X"])) {
    if(Rf_inherits(dataModel["X"], "dgCMatrix")) {
      X_ = as<sp_mat>(dataModel["X"]) ;
    }
    stop("unknown class of X") ;
  } else {
      X_ = as<mat>(dataModel["X"]) ;
  }
  
  standardize(center, scale) ;
}

template <typename matrix>
void RegressionData<matrix>::standardize(const bool center, const bool scale) {
  
  // TODO: OBSERVATION WEIGHTS
  
  if (center) {
    X_bar_ = mean(X_, 0).t();
    y_bar_ = mean(y_) ;
    // do not center X to preserve the sparsity structure
  } else {
    X_bar_ = zeros(p_) ;
    y_bar_ = 0;
  }
  
  if (scale) {
    norm_X_ = sqrt(trans(sum(square(X_))) - n_ * square(X_bar_));
    for (uword i=0; i<p_; i++) {
      X_.col(i) /= norm_X_(i);
    }
    X_bar_ /= norm_X_ ;
  } else {
    norm_X_ = ones(p_);
  }
  norm_y_ = sqrt(sum(square(y_))) ;
  
  if (center) {
    XTy_ = X_.t() * (y_- y_bar_) ;
    for (uword i=0; i<p_; i++) {
      XTy_(i) -=  sum(y_-y_bar_) * X_bar_(i);
    }
  } else {
    XTy_ = X_.t() * y_ ;
  }
  
  centered_ = center ;
  scaled_   = scale  ;
  
}

#endif
