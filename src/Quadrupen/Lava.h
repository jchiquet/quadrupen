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
  
  Lava(const RegressionData<matrix>&, const mat&, const bool&, const List&, const List&);

  void post_treatment(const RegressionData<matrix>& data,
                      const mat& U, const mat& V, const vec& D, const mat& C_inv) ;

  const sp_mat sparse_coefficients() const { 
    return sp_mat(join_cols(this->iA_, this->jA_), this->nzeros_, this->lambdas_.size(), this->data_.p_) ; 
  }

  const vector<vec>& dense_coefficients() const { 
    return dense_coef_ ; 
  }
  
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
  const bool& intercept,
  const List& regParam,
  const List& control) :
  ElasticNet<matrix>::ElasticNet(data, intercept, regParam, control), 
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
  
  sp_mat sparse_coef = this->coefficients() ;
  
  mat M = C_inv * V * diagmat(D/(square(D) + 1)) ;
  mat UTy = U.t() * data.y_ ;
  mat DVT = diagmat(D) * V.t();

  for (uword i=0;i< sparse_coef.n_rows ;i++){
    // Renormalize the sparse coefficient
    sparse_coef.row(i) /= data.norm_X_.t() ;
    // Compute the dense coefficient
    vec beta_ = M * (UTy - DVT * sparse_coef.row(i).t()) / data.norm_X_ ;
    dense_coef_.push_back(beta_) ;
    // Compute the constant / intercept term
    this->const_.push_back(data.y_bar_ - dot(sparse_coef.row(i).t() - beta_, data.X_bar_));
  }
  sparse_coef_ = sparse_coef ;
  
}

#endif
