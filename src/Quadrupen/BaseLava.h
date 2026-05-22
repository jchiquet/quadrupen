/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

// ====================================================
// Common base for Lava and GroupLava.
// Template on ConcreteBase (SparseRegularizer or GroupSparseRegularizer).
// Encapsulates:
//   - Lava preprocessing (SVD / K12 decomposition)
//   - post_treatment() (dense component + sparse refit)
//   - get_df() shared implementation

#pragma once

#include "RegressionData.h"

using arma::vec;
using arma::mat;
using arma::sp_mat;
using arma::uvec;
using arma::uword;
using arma::ones;
using Rcpp::List;

template <typename matrix, typename ConcreteBase>
class BaseLava : public ConcreteBase {

protected:

  // ── Preprocessing ─────────────────────────────────────────────────────────────
  struct LavaData {
    RegressionData<matrix> scaled ;
    mat Proj   ;
    mat B_proj ;
  } ;

  static LavaData lava_preprocess(RegressionData<matrix> orig, double gamma) {
    orig.scale_struct(gamma) ;
    mat Xc = orig.X_ ; Xc.each_row() -= orig.X_bar_.t() ;
    mat C_inv = solve(trimatu(chol(orig.S_.as_dense())), eye(orig.p_, orig.p_)) ;
    mat U, V ; vec D ;
    svd_econ(U, D, V, Xc * C_inv) ;
    vec D2     = square(D) ;
    mat Proj   = U * diagmat(D2 / (D2 + 1))       * U.t() ;
    mat K12    = U * diagmat(1 / sqrt(D2 + 1))     * U.t() ;
    mat B_proj = C_inv * V * diagmat(D / (D2 + 1)) * U.t() ;
    return {
      RegressionData<matrix>(K12 * Xc, K12 * (orig.y_ - orig.y_bar_),
                             sp_mat(orig.p_, orig.p_), ones(orig.p_), false, false),
      std::move(Proj),
      std::move(B_proj)
    } ;
  }

public:

  // ── Members set at construction ──────────────────────────────────────────────
  RegressionData<matrix> orig_data_ ;
  mat Proj_   ;
  mat B_proj_ ;

  // ── Members set by post_treatment ────────────────────────────────────────────
  sp_mat sparse_coef_ ;
  mat    b_           ;

  // Variadic constructor: forwards ConcreteBase-specific args after LavaData.
  // Lava:      Base(orig, t, regParam, control)
  // GroupLava: Base(orig, t, group_ind, regParam, control)
  template <typename... BaseArgs>
  BaseLava(const RegressionData<matrix>& orig_data, LavaData t, BaseArgs&&... base_args)
    : ConcreteBase(t.scaled, std::forward<BaseArgs>(base_args)...),
      orig_data_(orig_data),
      Proj_(std::move(t.Proj)),
      B_proj_(std::move(t.B_proj)) {}

  // ── get_df: base df from ConcreteBase + Lava projection correction ───────────
  double get_df() override {
    double df = ConcreteBase::get_df() ;
    this->set_.inverse_Gram() ;
    mat K = diagmat(ones(this->data_.n_)) -
      this->data_.X_.cols(this->set_.A_) * this->set_.XATXAinv_ * this->data_.X_.cols(this->set_.A_).t() ;
    df -= trace(K * Proj_) ;
    return df ;
  }

  // ── post_treatment: compute dense component b, then refit sparse ──────────────
  void post_treatment() {
    mat Xs = orig_data_.X_ ; Xs.each_row() -= orig_data_.X_bar_.t() ;
    sp_mat beta = this->coefficients() ;
    b_ = B_proj_ * (
      (orig_data_.y_ - orig_data_.y_bar_) * ones(1, this->lambdas_.size()) -
      Xs * beta
    ) ;
    this->coef_ = diagmat(1 / orig_data_.norm_X_) * (b_ + beta.as_dense()) ;
    this->intercept_.clear() ;
    this->intercept_debiased_.clear() ;
    this->debiased_.clear() ;
    for (uword i = 0; i < this->lambdas_.size(); i++) {
      this->intercept_.push_back(orig_data_.y_bar_ - dot(beta.col(i) + b_.col(i), orig_data_.X_bar_)) ;
      uvec A = this->active_[i] ;
      mat XsA = Xs.cols(A) ; mat XsAT = XsA.t() ;
      vec w = (orig_data_.y_ - this->data_.y_bar_) - Xs * b_.col(i) ;
      vec beta_debiased = solve(XsAT * XsA, XsAT * w) ;
      this->debiased_.insert(this->debiased_.end(), beta_debiased.begin(), beta_debiased.end()) ;
      this->intercept_debiased_.push_back(
        this->data_.y_bar_ - dot(beta_debiased, this->data_.X_bar_(A)) - dot(b_.col(i), this->data_.X_bar_)) ;
      beta.col(i) /= orig_data_.norm_X_ ;
    }
    sparse_coef_ = beta ;
  }

} ;
