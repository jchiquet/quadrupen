/*
 * Author: Julien CHIQUET, INRAE
 *         julien.chiquet@inrae.fr
 *         Statistique et Génome
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_REGRESSION_DATA_H
#define _quadrupen_REGRESSION_DATA_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Use template to handle dense or sparse encoding (mat/sp_mat in armadillo)
template <typename matrix>
class RegressionData {

public:

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
  
  RegressionData() ;
  RegressionData(const Environment& dataModel, const bool& center, const bool& scale);
  RegressionData(const matrix X,  const vec y, const sp_mat S, const vec weights, const bool center, const bool scale) ;

  void standardize();

  void scale_struct(vec weights) ;

  void scale_regressors(vec weights) ;

};

template <typename matrix>
RegressionData<matrix>::RegressionData(
  const Environment& dataModel,
  const bool& center, const bool& scale):
  y_       (as<vec>    (dataModel["y"])) , // response vector
  S_       (as<sp_mat> (dataModel["S"])) , // structuring matrix
  weights_ (as<vec>    (dataModel["wy"])), // observation weights
  centered_(center), scaled_(scale)        // standardization options
{
  X_ = as<matrix>(dataModel["X"]) ;
  n_ = X_.n_rows ; p_ = X_.n_cols ;
  standardize() ;
}

// Constructor from arma objects
template <typename matrix>
RegressionData<matrix>::RegressionData(
  const matrix X, 
  const vec y, 
  const sp_mat S, 
  const vec weights,
  const bool center, const bool scale):
  X_       (X)    , // response vector
  y_       (y) , // response vector
  S_       (S) , // structuring matrix
  weights_ (weights), // observation weights
  centered_(center), scaled_(scale) // standardization options
{ n_ = X_.n_rows ; p_ = X_.n_cols ; }

template <typename matrix>
void RegressionData<matrix>::scale_struct(vec weights) {
  sp_mat diag_S = spdiags(weights, ivec({0}), p_, p_) ;
  S_ = diag_S * S_ * diag_S  ; 
};

template <>
inline void RegressionData<sp_mat>::scale_regressors(vec weights) {
  
  if (any(weights != 1)) {
    for (uword i=0; i<n_; i++) {
      X_.row(i) /= trans(weights) ;
    }
    X_bar_ /= weights ;
  }
  XTy_ = X_.t() * (y_-y_bar_) - sum(y_-y_bar_) * X_bar_ ;

};

template <>
inline void RegressionData<mat>::scale_regressors(vec weights) {
  
  if (any(weights != 1)) {
    X_.each_row() /= trans(weights) ;
    X_bar_ /= weights ;
  }
  XTy_ = X_.t() * (y_-y_bar_) - sum(y_-y_bar_) * X_bar_ ;

};

template <typename matrix>
void RegressionData<matrix>::standardize() {
  
  // TODO: OBSERVATION WEIGHTS
  
  if (centered_) {
    X_bar_ = mean(X_, 0).t();
    y_bar_ = mean(y_) ;
    // do not center X to preserve the sparsity structure
  } else {
    X_bar_ = zeros(p_) ;
    y_bar_ = 0;
  }
  
  if (scaled_) {
    norm_X_ = sqrt(trans(sum(square(X_))) - n_ * square(X_bar_));
    X_ = X_ * diagmat(1/norm_X_) ;
    X_bar_ /= norm_X_ ;
  } else {
    norm_X_ = ones(p_);
  }
  norm_y_ = sqrt(sum(square(y_))) ;
  
  XTy_ = X_.t() * (y_-y_bar_) - sum(y_-y_bar_) * X_bar_ ;
}

#endif
