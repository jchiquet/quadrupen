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
  
  Lava(const RegressionData<matrix>&, const mat&, const List&, const List&);

  void post_treatment(const RegressionData<matrix>& data,
                      const mat& U, const mat& V, const vec& D, const mat& C_inv) ;

  // Specific to Elastic-Net regularization
  sp_mat sparse_coef_     ; // matrix of dense coefficients
  vector<vec> dense_coef_ ; // matrix of dense coefficients
  vector<double> const_   ; // matrix of dense coefficients
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
  
  double df = this->set_.size() ;
  mat B ;
  if (this->set_.use_chol_) {
    B = this->set_.Rinv_ * this->set_.Rinv_.t();
  } else {
    B = inv_sympd(this->set_.XATXA_);
  }
  mat K = diagmat(ones(this->data_.n_)) - 
    this->data_.X_.cols(this->set_.A_) * B * this->data_.X_.cols(this->set_.A_).t() ;
  
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
    this->const_.push_back(data.y_bar_ - dot(beta.row(i).t() - b, data.X_bar_));
  }
  sparse_coef_ = beta ;
  
}

#endif
