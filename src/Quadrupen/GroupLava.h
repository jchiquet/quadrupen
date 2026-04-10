/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _Group_Lava_H
#define _Group_Lava_H

#include "GroupElasticNet.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, MixedNorm norm>
class GroupLava : public GroupElasticNet<matrix,norm> {
  
public:
  
  using Regularizer<matrix>::coef_       ;
  using Regularizer<matrix>::intercept_  ;
  using Regularizer<matrix>::data_       ;
  using GroupElasticNet<matrix,norm>::set_     ;
  using GroupElasticNet<matrix,norm>::active_  ;
  using GroupElasticNet<matrix,norm>::debiased_ ;
  using GroupElasticNet<matrix,norm>::intercept_debiased_ ;
  using GroupElasticNet<matrix,norm>::lambda_factor_ ;
  
  GroupLava(const RegressionData<matrix>&, const mat&, const uvec&, const List&, const List&);
  
  void post_treatment(const RegressionData<matrix>& data, const mat& b) ;
  
  // Specific to LAVA
  sp_mat sparse_coef_ ; // matrix of sparse coefficients
  mat Proj_           ; // Lava projector matrix
  
  // Compute degrees of freedom for the current estimate
  double get_df() ;
  
};

template <typename matrix, MixedNorm norm>
GroupLava<matrix,norm>::GroupLava(
  const RegressionData<matrix>& data,
  const mat& Proj,
  const uvec& group_ind, 
  const List& regParam,
  const List& control) :
  GroupElasticNet<matrix,norm>::GroupElasticNet(data, group_ind, regParam, control), 
  Proj_(Proj) {}

template <typename matrix, MixedNorm norm>
double GroupLava<matrix,norm>::get_df() {

  double df = GroupElasticNet<matrix,norm>::get_df() ;
  mat K = diagmat(ones(data_.n_)) - 
    data_.X_.cols(set_.A_) * set_.XATXAinv_ * data_.X_.cols(set_.A_).t() ;
  
  df -= trace(K * Proj_);
  
  return(df);
}

template <typename matrix, MixedNorm norm>
void GroupLava<matrix,norm>::post_treatment(const RegressionData<matrix>& data, const mat& b) {
  
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
    debiased_ = join_cols(debiased_, beta_debiased/(data_.norm_X_(A)));
    // debiased_ = join_cols(debiased_, beta_debiased/(data_.norm_X_(A) % lambda_factor_(A)));
    intercept_debiased_.push_back(
      data_.y_bar_ -
        dot(beta_debiased, data_.X_bar_(A)) - dot(b.col(i), data_.X_bar_)) ;
    
    // Scale the sparse coefficients to original
    beta.col(i) /= data.norm_X_ ;
  }
  
  sparse_coef_ = beta ;
  
}

#endif
