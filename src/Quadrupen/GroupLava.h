/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "GroupSparseRegularizer.h"

using arma::vec;
using arma::mat;
using arma::sp_mat;
using arma::uvec;
using arma::uword;
using arma::ones;
using Rcpp::List;
using std::vector;

template <typename matrix, GroupSparseNorm norm>
class GroupLava : public GroupSparseRegularizer<matrix,norm> {

public:

  using Regularizer<matrix>::coef_       ;
  using Regularizer<matrix>::intercept_  ;
  using Regularizer<matrix>::data_       ;
  using GroupSparseRegularizer<matrix,norm>::set_     ;
  using GroupSparseRegularizer<matrix,norm>::active_  ;
  using GroupSparseRegularizer<matrix,norm>::debiased_ ;
  using GroupSparseRegularizer<matrix,norm>::intercept_debiased_ ;
  using GroupSparseRegularizer<matrix,norm>::lambda_factor_ ;

  GroupLava(const RegressionData<matrix>&, const mat&, const uvec&, const List&, const List&);

  void post_treatment(const RegressionData<matrix>& data, const mat& b) ;

  // Specific to LAVA
  sp_mat sparse_coef_ ; // matrix of sparse coefficients
  mat Proj_           ; // Lava projector matrix

  // Compute degrees of freedom for the current estimate
  double get_df() override ;

};

template <typename matrix, GroupSparseNorm norm>
GroupLava<matrix,norm>::GroupLava(
    const RegressionData<matrix>& data,
    const mat& Proj,
    const uvec& group_ind,
    const List& regParam,
    const List& control) :
  GroupSparseRegularizer<matrix,norm>::GroupSparseRegularizer(data, group_ind, regParam, control),
  Proj_(Proj) {}

template <typename matrix, GroupSparseNorm norm>
double GroupLava<matrix,norm>::get_df() {

  double df = GroupSparseRegularizer<matrix,norm>::get_df() ;
  // Materialise XATXAinv_ here: store_path_step() now uses solve_Gram() and skips this
  set_.inverse_Gram() ;
  mat K = diagmat(ones(data_.n_)) -
    data_.X_.cols(set_.A_) * set_.XATXAinv_ * data_.X_.cols(set_.A_).t() ;

  df -= trace(K * Proj_);

  return(df);
}

template <typename matrix, GroupSparseNorm norm>
void GroupLava<matrix,norm>::post_treatment(const RegressionData<matrix>& orig_data, const mat& b) {
  this->post_treatment_impl(orig_data, b, sparse_coef_,
    [this](const uvec& A) { return data_.norm_X_(A) ; }) ;
}
