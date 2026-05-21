/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

// ====================================================
// Intermediate base class shared by SparseRegularizer and GroupSparseRegularizer.
// Holds the common state (nzeros_, debiased_, active_, …) and implements
// solution_path() via three pure-virtual hooks:
//   - run_solver(lambda)  : calls working_set on the concrete solver/active-set
//   - store_path_step()   : stores nzeros, intercept, debiased, active for this lambda
//   - solver_diagnostics(): returns (inner_iter, J_vec, D_vec) for the output List

#pragma once

#include "Regularizer.h"
#include "ActiveSet.h"

using arma::vec;
using arma::uvec;
using arma::umat;
using arma::uword;
using arma::sp_mat;
using arma::wall_clock;
using Rcpp::List;
using Rcpp::Named;
using std::vector;

template <typename matrix>
class SparsifyingRegularizer : public Regularizer<matrix> {

public:

  using Regularizer<matrix>::Regularizer ;       // inherit constructors
  using Regularizer<matrix>::intercept_     ;
  using Regularizer<matrix>::lambdas_       ;
  using Regularizer<matrix>::gamma_         ;
  using Regularizer<matrix>::data_          ;
  using Regularizer<matrix>::df_            ;
  using Regularizer<matrix>::beta_          ;
  using Regularizer<matrix>::grad_          ;
  using Regularizer<matrix>::lambda_factor_ ;
  using Regularizer<matrix>::get_lambda_seq ;
  using Regularizer<matrix>::build_sp_locations ;

  // ── Common state (all lambdas) ──────────────────────────────────────────────
  std::vector<double> nzeros_             ; // scaled nonzero coefficients
  std::vector<double> debiased_           ; // scaled debiased coefficients
  vector<uvec>        active_             ; // active index sets
  vector<double>      intercept_debiased_ ; // debiased intercepts
  vec                 beta_debiased_      ; // current debiased beta

  // ── Accessors ───────────────────────────────────────────────────────────────
  const vector<double>& intercept_debiased() const { return intercept_debiased_ ; }

  const sp_mat coefficients() const {
    return sp_mat(build_sp_locations(active_),
                  vec(nzeros_), data_.p_, active_.size(), true, false) ;
  }

  const sp_mat debiased_coefficients() const {
    return sp_mat(build_sp_locations(active_),
                  vec(debiased_), data_.p_, active_.size(), true, false) ;
  }

  const sp_mat active_var() const {
    umat locs = build_sp_locations(active_) ;
    return sp_mat(locs, arma::ones<vec>(locs.n_cols),
                  data_.p_, active_.size(), true, false) ;
  }

  // ── Pure-virtual hooks implemented by each subclass ────────────────

  // ── Degrees of freedom ──────────────────────────────────────────────────────
  virtual double get_df() = 0 ;

  // ── Active set accessor (implemented by each concrete subclass) ─────────────
  virtual ActiveSet<matrix>& current_set() = 0 ;

  // Call working_set for the given lambda; return (status, gap, outer_iter).
  struct StepResult { uword status ; double gap ; uword iter ; } ;
  virtual StepResult run_solver(double lambda) = 0 ;

  // Return the solver's inner diagnostic vectors for the output List.
  struct Diagnostics { const vector<uword>& inner_iter ; const vector<double>& J_vec ; const vector<double>& D_vec ; } ;
  virtual Diagnostics solver_diagnostics() const = 0 ;

  // ── Shared solution_path ────────────────────────────────────────────────────
  List solution_path(const List& /*control*/) {

    vector<double> gap, timing ;
    vector<uword>  status, iter ;
    wall_clock timer ; timer.tic() ;

    for (auto lambda_ : lambdas_) {

      StepResult step = run_solver(lambda_) ;
      status.push_back(step.status) ;
      gap.push_back(step.gap) ;
      iter.push_back(step.iter) ;

      if (status.back() >= 2) {
        break ;
      } else {
        auto& set = current_set() ;
        vec nz = beta_ / data_.norm_X_(set.A_) ;
        nzeros_.insert(nzeros_.end(), nz.begin(), nz.end()) ;
        intercept_.push_back(data_.y_bar_ - dot(beta_, data_.X_bar_(set.A_))) ;
        beta_debiased_ = set.solve_Gram(data_.XTy_(set.A_)) ;
        vec db = beta_debiased_ / data_.norm_X_(set.A_) ;
        debiased_.insert(debiased_.end(), db.begin(), db.end()) ;
        intercept_debiased_.push_back(data_.y_bar_ - dot(beta_debiased_, data_.X_bar_(set.A_))) ;
        active_.push_back(set.A_) ;
       df_.push_back(get_df()) ;
      }

      timing.push_back(timer.toc()) ;
    }
    lambdas_.resize(df_.size()) ;

    Diagnostics diag = solver_diagnostics() ;
    return List::create(
      Named("it_active")      = iter,
      Named("it_optim")       = diag.inner_iter,
      Named("max_grd")        = gap,
      Named("gap_hat")        = diag.J_vec,
      Named("delta_hat")      = diag.D_vec,
      Named("convergence")    = status,
      Named("pensteps_timer") = timing
    ) ;
  }

} ;
