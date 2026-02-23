/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _Lava_H
#define _Lava_H

#include "ElasticNet.h"
#include "OptimizerL1.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class Lava : 
  public ElasticNet<matrix> {

public:

  using GenericRegularizer<matrix,Norm::L1>::intercept_  ;
  using GenericRegularizer<matrix,Norm::L1>::set_        ;
  using GenericRegularizer<matrix,Norm::L1>::data_       ;
  using ElasticNet<matrix>::debiased_ ;
  using ElasticNet<matrix>::intercept_debiased_ ;
  using ElasticNet<matrix>::lambda_factor_ ;
  
  Lava(const RegressionData<matrix>&, const mat&, const List&, const List&);

  void post_treatment(const RegressionData<matrix>& data,
                      const mat& U, const mat& V, const vec& D, const mat& C_inv) ;

  // Specific to Elastic-Net regularization
  sp_mat sparse_coef_     ; // matrix of dense coefficients
  vector<vec> dense_coef_ ; // matrix of dense coefficients
  mat Proj_               ; // Lava projector matrix

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
  
  double df = set_.size() ;
  mat K = diagmat(ones(data_.n_)) - 
    data_.X_.cols(set_.A_) * set_.XATXAinv_ * data_.X_.cols(set_.A_).t() ;
  
  df -= trace(K * Proj_);
  
  return(df);
}

template <typename matrix>
void Lava<matrix>::post_treatment(const RegressionData<matrix>& data,
                                  const mat& U, const mat& V, const vec& D, const mat& C_inv) {
  
  sp_mat beta = this->coefficients() ;
  
  mat M = C_inv * V * diagmat(D/(square(D) + 1)) ;
  mat UTy = U.t() * (data.y_ - data.y_bar_) ;
  mat DVT = diagmat(D) * V.t();
  debiased_.reset() ;
  mat X = data.X_ ;
  X.each_row() -= data.X_bar_.t() ;
  
  for (uword i=0;i< beta.n_rows ;i++){
    // Rescale sparse and dense coefficient to original scaling factors
    beta.row(i) /= data.norm_X_.t() ;
    vec b = M * (UTy - DVT * beta.row(i).t()) / data.norm_X_ ;
    dense_coef_.push_back(b) ;
    intercept_.push_back(data.y_bar_ - dot(beta.row(i).t() - b, data.X_bar_));
    
    // Refit the sparse coefficient to remove bias due to sparse shrinkage
    uvec A = find(beta.row(i)) ;
    vec w = data.y_ - data.X_ * b  - dot(data.X_bar_,b) * ones(data.n_);
    vec beta_debiased = solve(X.cols(A).t() * X.cols(A), X.cols(A).t() * w) ;
    debiased_ = join_cols(debiased_, beta_debiased/(data_.norm_X_(A) % lambda_factor_(A)));
    intercept_debiased_.push_back(data_.y_bar_ - dot(beta_debiased, data_.X_bar_(A))) ;
  }

  sparse_coef_ = beta ;

}

#endif
