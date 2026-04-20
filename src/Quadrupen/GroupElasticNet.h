/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _GroupLasso_H
#define _GroupLasso_H

#include "RegularizerSparse.h"
#include "OptimizerGroup.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, GroupNorm norm>
class GroupElasticNet : 
  public GroupSparseRegularizer<matrix,norm>{
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
    using SparseRegularizer<matrix>::nzeros_   ;
    using SparseRegularizer<matrix>::debiased_ ;
    using SparseRegularizer<matrix>::beta_debiased_ ;
    using SparseRegularizer<matrix>::intercept_debiased_   ;
    using SparseRegularizer<matrix>::active_ ;
    using GroupSparseRegularizer<matrix,norm>::penalty_   ;
    using GroupSparseRegularizer<matrix,norm>::set_  ;
    using GroupSparseRegularizer<matrix,norm>::get_lambda_max ;

    GroupElasticNet(const RegressionData<matrix>&, const uvec&, const List&, const List&);
  
    List solution_path(const List&);
    
    // Specific to Group Lasso regularization
    GroupOptimizer<matrix,norm> solver_ ; // Solvers for Group L1 penalty

    // Compute degrees of freedom for the current estimate
    double get_df() ;
    
};

template <typename matrix, GroupNorm norm>
GroupElasticNet<matrix,norm>::GroupElasticNet(
  const RegressionData<matrix>& data, const uvec& group_ind, const List& regParam, const List& control) :
  GroupSparseRegularizer<matrix,norm>::GroupSparseRegularizer(data, regParam) {

    // Scale the structuring matrix according to the amount of l2 penalty 
    data_.scale_struct(gamma_) ;
    
    // Initialize the active set, beta_ and gradient with starting coefficient
    vec beta0 = control["beta0"] ;
    // uvec A0 = find(beta0) ;
    // if (A0.is_empty()) {
    set_ = ActiveSetGroup(data_, group_ind, as<bool>(control["usechol"])) ;
    grad_ = - data_.XTy_ ;
    // } else {
    //   set_  = ActiveSet(data_, A0, as<bool>(control["usechol"])) ;
    //   beta_ = beta0(A0) ;
    //   grad_ = - data_.XTy_ + set_.XTXA_ * beta_  ;
    // }

    // set up the penalty
    penalty_ = GroupPenalty<norm>(as<double>(regParam["alpha"])) ;

    get_lambda_seq(get_lambda_max(), regParam) ;
    
    // Set up the optimizer
    solver_ = GroupOptimizer<matrix,norm>(penalty_, control);

  }

template <typename matrix, GroupNorm norm>
double GroupElasticNet<matrix,norm>::get_df() {

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

template <typename matrix, GroupNorm norm>
List GroupElasticNet<matrix,norm>::solution_path(const List& control) {
  
  vector<double> gap, timing ; // timings and optimality measures
  vector<uword> status, iter ; // convergence and # of inner/outer iterates
  
  // LAMBDA LOOP
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda_ : lambdas_) {

    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    status.push_back(
      solver_.solve(beta_, grad_, lambda_, lambda_factor_, gamma_, data_, set_)
    ) ;
    gap.push_back(solver_.gap_) ;
    iter.push_back(solver_.iter_) ;

    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      // store current coefficients
      nzeros_ = join_cols(nzeros_, beta_/data_.norm_X_(set_.A_));
      intercept_.push_back(data_.y_bar_ - dot(beta_, data_.X_bar_(set_.A_))); // X_bar is scaled
      // compute and store debiased coefficients
      set_.inverse_Gram() ;
      beta_debiased_ = set_.XATXAinv_ * (data_.XTy_(set_.A_) - data_.X_bar_(set_.A_) * accu(data_.y_)) ;
      debiased_ = join_cols(debiased_, beta_debiased_/data_.norm_X_(set_.A_));
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

#endif

