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
  mat B ;
  if (set_.use_chol_) {
    B = set_.Rinv_ * set_.Rinv_.t();
  } else {
    B = inv_sympd(set_.XATXA_);
  }
  mat K = diagmat(ones(data_.n_)) - 
    data_.X_.cols(set_.A_) * B * data_.X_.cols(set_.A_).t() ;
  
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

  for (uword i=0;i< beta.n_rows ;i++){
    beta.row(i) /= data.norm_X_.t() ;
    vec b = M * (UTy - DVT * beta.row(i).t()) / data.norm_X_ ;
    dense_coef_.push_back(b) ;
    intercept_.push_back(data.y_bar_ - dot(beta.row(i).t() - b, data.X_bar_));
  }
  sparse_coef_ = beta ;
  
}

#endif
