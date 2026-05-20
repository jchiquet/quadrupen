/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

// ====================================================
// Group-Sparse Regularizers       

#pragma once

#include "Regularizer.h"
#include "OptimizerGroup.h"

using arma::vec;
using arma::uvec;
using arma::uword;
using arma::sp_mat;
using arma::urowvec;
using arma::wall_clock;
using arma::join_cols;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::as;
using std::vector;

template <typename matrix, GroupSparseNorm norm>
class GroupSparseRegularizer : 
  public Regularizer<matrix> {
  public:

    using Regularizer<matrix>::intercept_ ;
    using Regularizer<matrix>::lambdas_   ;
    using Regularizer<matrix>::gamma_   ;
    using Regularizer<matrix>::data_ ;
    using Regularizer<matrix>::df_   ;
    using Regularizer<matrix>::beta_   ;
    using Regularizer<matrix>::grad_   ;
    using Regularizer<matrix>::lambda_factor_ ;
    using Regularizer<matrix>::get_lambda_seq ;
    
    std::vector<double> nzeros_  ; // contains non-zero values of all betas (for all lambda values)
    std::vector<double> debiased_ ; // contains debiased non-zero values of all betas (for all lambda values)
    vector<uvec> active_ ; // successively activated variable (for all lambda values)
    vector<double >intercept_debiased_ ; // debiased vector of intercept values (for all lambda values)
    vec beta_debiased_ ; // vector of current active beta debiased (for the current lambda)
    ActiveSetGroup<matrix> set_ ; // Active set of variable and data
    GroupPenalty<norm> penalty_ ; // main penalty object 
    
    double get_lambda_max() 
    {
      return(
        penalty_.lambda_max(data_.XTy_, set_.grp_sizes_, lambda_factor_)
      );
    }

    GroupSparseRegularizer(const RegressionData<matrix>&, const uvec&, const List&, const List&);
  
    List solution_path(const List&);
    
    // Specific to Group Lasso regularization
    GroupOptimizer<matrix,norm> solver_ ; // Solvers for Group Sparse penalty

    // Compute degrees of freedom for the current estimate
    double get_df() ;

    const sp_mat coefficients() const {
      vector<uword> rowA, colA ;
      uword current_col = 0;
      for (const auto& a : active_) {
        rowA.insert(rowA.end(), a.begin(), a.end()) ;
        colA.insert(colA.end(), a.n_elem, current_col) ;
        current_col++;
      }
      return sp_mat(join_cols(urowvec(rowA), urowvec(colA)),
                    vec(nzeros_), data_.p_, active_.size(), true, false) ;
    }

    const sp_mat debiased_coefficients() const {
      vector<uword> rowA, colA ;
      uword current_col = 0;
      for (const auto& a : active_) {
        rowA.insert(rowA.end(), a.begin(), a.end()) ;
        colA.insert(colA.end(), a.n_elem, current_col) ;
        current_col++;
      }
      return sp_mat(join_cols(urowvec(rowA), urowvec(colA)),
                    vec(debiased_), data_.p_, active_.size(), true, false) ;
    }
    
    const sp_mat active_var() const { 
      vector<uword> rowA, colA ;
      uword current_col = 0;
      for (const auto& a : active_) {
        rowA.insert(rowA.end(), a.begin(), a.end()) ;
        colA.insert(colA.end(), a.n_elem, current_col) ;
        current_col++;
      }
      return sp_mat(join_cols(urowvec(rowA), urowvec(colA)),
                    vec(rowA.size(), arma::fill::ones), data_.p_, active_.size(), true, false) ;
    }
    
    const vector<double>& intercept_debiased() const { return intercept_debiased_ ; }
    
};

template <typename matrix, GroupSparseNorm norm>
GroupSparseRegularizer<matrix,norm>::GroupSparseRegularizer(
  const RegressionData<matrix>& data, const uvec& group_ind, const List& regParam, const List& control) :
  Regularizer<matrix>::Regularizer(data, regParam) {

    // Scale the structuring matrix according to the amount of l2 penalty 
    data_.scale_struct(gamma_) ;
    
    // Initialize the active set and the gradient
    grad_ = - data_.XTy_ ;
    set_ = ActiveSetGroup(data_, group_ind, as<bool>(control["factmat"])) ;
    
    // Set up the penalty
    penalty_ = GroupPenalty<norm>(as<double>(regParam["alpha"])) ;
    get_lambda_seq(penalty_.lambda_max(data_.XTy_, set_.grp_sizes_, lambda_factor_), regParam) ;

    // Set up the optimizer
    solver_ = GroupOptimizer<matrix,norm>(penalty_, control);

  }

template <typename matrix, GroupSparseNorm norm>
double GroupSparseRegularizer<matrix,norm>::get_df() {

  double df = data_.centered_ ;

  if (set_.size_grp() > 0) {
    // approximate degrees of freedom
    vec active_grp_norm = penalty_.elt_norm(beta_, set_.grp_sizes_(set_.G_), ones(set_.size_grp())) ;
    vec active_grp_norm_ols = penalty_.elt_norm(beta_debiased_, set_.grp_sizes_(set_.G_), ones(set_.size_grp())) / (1 + gamma_) ;
    
    df = df + 
      accu(1 + (active_grp_norm / active_grp_norm_ols) % (set_.grp_sizes_(set_.G_) - 1)) ;
  }  

  return(df);
}

template <typename matrix, GroupSparseNorm norm>
List GroupSparseRegularizer<matrix,norm>::solution_path(const List& control) {
  
  vector<double> gap, timing ; // timings and optimality measures
  vector<uword> status, iter ; // convergence and # of inner/outer iterates
  // LAMBDA LOOP
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda_ : lambdas_) {

    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    status.push_back(
      solver_.working_set(beta_, grad_, lambda_, lambda_factor_, gamma_, data_, set_)
    ) ;
    gap.push_back(solver_.gap_) ;
    iter.push_back(solver_.iter_) ;

    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      // store current coefficients
      vec nz = beta_ / data_.norm_X_(set_.A_) ;
      nzeros_.insert(nzeros_.end(), nz.begin(), nz.end()) ;
      intercept_.push_back(data_.y_bar_ - dot(beta_, data_.X_bar_(set_.A_))); // X_bar is scaled
      // compute and store debiased coefficients
      set_.inverse_Gram() ;
      // XTy_(A) = X_A'(y - y_bar) already incorporates centering; no further correction needed
      beta_debiased_ = set_.XATXAinv_ * data_.XTy_(set_.A_) ;
      vec db = beta_debiased_ / data_.norm_X_(set_.A_) ;
      debiased_.insert(debiased_.end(), db.begin(), db.end()) ;
      intercept_debiased_.push_back(data_.y_bar_ - dot(beta_debiased_, data_.X_bar_(set_.A_))) ; // X_bar is scaled
      // store degrees fo freedom and current active set
      df_.push_back(this->get_df()) ;
      active_.push_back(set_.A_) ;
    }

    timing.push_back(timer.toc()) ;
  } // END OF THE LOOP OVER LAMBDA
  lambdas_.resize(df_.size()) ;
  
  return(
    List::create(
      Named("it_active")      = iter,
      Named("it_optim")       = solver_.inner_iter_ ,
      Named("max_grd")        = gap,
      Named("gap_hat")        = solver_.J_vec_,
      Named("delta_hat")      = solver_.D_vec_,
      Named("convergence")    = status ,
      Named("pensteps_timer") = timing
    )
  );
}


