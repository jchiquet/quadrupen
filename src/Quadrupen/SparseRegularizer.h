/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */
#pragma once

#include "SparsifyingRegularizer.h"
#include "OptimizerSparse.h"

using arma::vec;
using arma::uvec;
using arma::uword;
using arma::sp_mat;
using arma::umat;
using Rcpp::List;
using Rcpp::as;
using std::vector;

template <typename matrix, SparseNorm norm>
class SparseRegularizer :
  public SparsifyingRegularizer<matrix> {
  public:

    using SparsifyingRegularizer<matrix>::intercept_ ;
    using SparsifyingRegularizer<matrix>::gamma_     ;
    using SparsifyingRegularizer<matrix>::data_      ;
    using SparsifyingRegularizer<matrix>::beta_      ;
    using SparsifyingRegularizer<matrix>::grad_      ;
    using SparsifyingRegularizer<matrix>::lambda_factor_ ;
    using SparsifyingRegularizer<matrix>::get_lambda_seq ;
    using SparsifyingRegularizer<matrix>::nzeros_    ;
    using SparsifyingRegularizer<matrix>::debiased_  ;
    using SparsifyingRegularizer<matrix>::active_    ;
    using SparsifyingRegularizer<matrix>::intercept_debiased_ ;
    using SparsifyingRegularizer<matrix>::beta_debiased_ ;
    using typename SparsifyingRegularizer<matrix>::StepResult ;
    using typename SparsifyingRegularizer<matrix>::Diagnostics ;

    SparsePenalty<norm>          penalty_ ;
    ActiveSet<matrix>            set_     ;
    SparseOptimizer<matrix,norm> solver_  ;

    double get_lambda_max() {
      return penalty_.lambda_max(data_.XTy_, lambda_factor_) ;
    }

    SparseRegularizer(const RegressionData<matrix>&, const List&, const List&) ;

    double get_df() override ;

    // ── SparsifyingRegularizer hooks ─────────────────────────────────────────

    ActiveSet<matrix>& current_set() override { return set_ ; }

    StepResult run_solver(double lambda) override {
      uword status = solver_.working_set(beta_, grad_, lambda, lambda_factor_, gamma_, data_, set_) ;
      return { status, solver_.gap_, solver_.iter_ } ;
    }

    Diagnostics solver_diagnostics() const override {
      return { solver_.inner_iter_, solver_.J_vec_, solver_.D_vec_ } ;
    }

} ;

template <typename matrix, SparseNorm norm>
SparseRegularizer<matrix,norm>::SparseRegularizer(
  const RegressionData<matrix>& data, const List& regParam, const List& control) :
  SparsifyingRegularizer<matrix>::SparsifyingRegularizer(data, regParam) {

    // Scale the structuring matrix according to the amount of l2 penalty
    data_.scale_struct(gamma_) ;

    // Initialize the active set, beta_ and the gradient
    vec beta0 = control["beta0"] ;
    uvec A0 = find(beta0) ;
    grad_ = - data_.XTy_ ;
    if (A0.is_empty()) {
      set_  = ActiveSet(data_, as<bool>(control["factmat"])) ;
    } else {
      set_  = ActiveSet(data_, A0, as<bool>(control["factmat"])) ;
      beta_ = beta0(A0) ;
      grad_ += set_.XTXA_ * beta_ ;
    }

    // Set the penalty
    penalty_ = SparsePenalty<norm>(as<double>(regParam["eta"])) ;
    get_lambda_seq(penalty_.lambda_max(data_.XTy_, lambda_factor_), regParam) ;

    // Set up the optimizer
    solver_ = SparseOptimizer<matrix,norm>(penalty_, control) ;
  }

template <typename matrix, SparseNorm norm>
double SparseRegularizer<matrix,norm>::get_df() {

  double df = set_.size() + data_.centered_ ;
  if (gamma_ > 0) {
    // Inversion of XATXAinv_ (only paid when gamma_ > 0)
    set_.inverse_Gram() ;
    // loop due to sparse encoding. should iterate over the n_zeros only...
    mat SAA(set_.size(), set_.size()) ;
    for (uword i=0;i<set_.size();i++){
      for (uword j=i;j<set_.size();j++){
        SAA(i,j) = data_.S_.at(set_.A_(i),set_.A_(j));
        SAA(j,i) = SAA(i,j);
      }
    }
    // note that XATXAinv_ is in fact (XATXA + lambda S)^-1
    // trace(A*B) = sum(A .* B') avoids materialising the k×k product: O(k²) vs O(k³)
    df -= accu(SAA % set_.XATXAinv_);
  }

  return df ;
}
