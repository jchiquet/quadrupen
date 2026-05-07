/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "ElasticNet.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class Lava : public ElasticNet<matrix> {

public:

  using Regularizer<matrix>::coef_       ;
  using Regularizer<matrix>::intercept_  ;
  using Regularizer<matrix>::data_       ;
  using ElasticNet<matrix>::set_      ;
  using ElasticNet<matrix>::active_  ;
  using ElasticNet<matrix>::debiased_ ;
  using ElasticNet<matrix>::intercept_debiased_ ;
  using ElasticNet<matrix>::lambda_factor_ ;
  
  Lava(const RegressionData<matrix>&, const mat&, const List&, const List&);

  void post_treatment(const RegressionData<matrix>& data, const mat& b) ;

  // Specific to LAVA
  sp_mat sparse_coef_ ; // matrix of sparse coefficients
  mat Proj_           ; // Lava projector matrix

  // Compute degrees of freedom for the current estimate
  double get_df() ;
  
};

template <typename matrix>
Lava<matrix>::Lava(
  const RegressionData<matrix>& data,
  const mat& Proj,
  const List& regParam,
  const List& control) :
  ElasticNet<matrix>::ElasticNet(data, regParam, control), 
  Proj_(Proj) {}

template <typename matrix>
double Lava<matrix>::get_df() {
  
  double df = set_.size() + data_.centered_ ;
  mat K = diagmat(ones(data_.n_)) - 
    data_.X_.cols(set_.A_) * set_.XATXAinv_ * data_.X_.cols(set_.A_).t() ;
  
  df -= trace(K * Proj_);
  
  return(df);
}

template <typename matrix>
void Lava<matrix>::post_treatment(const RegressionData<matrix>& data, const mat& b) {
  
  sp_mat beta = this->coefficients() ;

  mat Xs = data.X_ ; Xs.each_row() -= data.X_bar_.t() ;
  
  // Scale coefficients (theta = beta + b) to original
  coef_ = diagmat(1/data.norm_X_) * (b + beta.as_dense()) ; 
  
  debiased_.reset() ;
  intercept_.clear() ;
  intercept_debiased_.clear() ;
  
  for (uword i=0; i< this->lambdas_.size() ;i++){
    
    intercept_.push_back(data.y_bar_ - dot(beta.col(i) + b.col(i), data.X_bar_));

    // Refit the sparse coefficient to remove bias due to sparse shrinkage
    uvec A = active_[i] ;
    vec w = (data.y_ - data_.y_bar_) - Xs * b.col(i) ;
    vec beta_debiased = solve(Xs.cols(A).t() * Xs.cols(A), Xs.cols(A).t() * w) ;
    debiased_ = join_cols(debiased_, beta_debiased/(data_.norm_X_(A) % lambda_factor_(A)));
    intercept_debiased_.push_back(
      data_.y_bar_ -
        dot(beta_debiased, data_.X_bar_(A)) - dot(b.col(i), data_.X_bar_)) ;
    
    // Scale the sparse coefficients to original
    beta.col(i) /= data.norm_X_ ;
  }

  sparse_coef_ = beta ;

}

