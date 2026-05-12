/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */
#pragma once

#include "Regularizer.h"
#include "OptimizerSparse.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, SparseNorm norm>
class SparseRegularizer : 
  public Regularizer<matrix>{
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
    
    vec   nzeros_ ; // contains non-zero values of all betas (for all lambda values)
    vec debiased_ ; // contains debiased non-zero values of all betas (for all lambda values)
    vector<uvec> active_ ; // successively activated variable (for all lambda values)
    vector<double >intercept_debiased_ ; // debiased vector of intercept values (for all lambda values)
    vec beta_debiased_ ; // vector of current active beta debiased (for the current lambda)
    SparsePenalty<norm> penalty_ ; // main penalty object 
    ActiveSet<matrix> set_       ; // Active set of variable and data
    
    double get_lambda_max() 
    {
      return(penalty_.lambda_max(data_.XTy_, lambda_factor_));
    }
    
    SparseRegularizer(const RegressionData<matrix>&, const List&, const List&);
    
    List solution_path(const List&);
    
    // Specific to Elastic-Net regularization
    SparseOptimizer<matrix, norm> solver_ ; // Solvers for L1 penalty
    
    void optimality_violation(double lambda_, uword type) ;
    
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
                    nzeros_, data_.p_, active_.size(), true, false) ;
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
                    debiased_, data_.p_, active_.size(), true, false) ;
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
                    vec(rowA.size(), fill::ones), data_.p_, active_.size(), true, false) ;
    }
    
    const vector<double>& intercept_debiased() const { return intercept_debiased_ ; }

};

template <typename matrix, SparseNorm norm>
SparseRegularizer<matrix,norm>::SparseRegularizer(
  const RegressionData<matrix>& data, const List& regParam, const List& control) :
  Regularizer<matrix>::Regularizer(data, regParam) {

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
      grad_ += set_.XTXA_ * beta_  ;
    }
    
    // Set the penalty to l1
    penalty_ = SparsePenalty<norm>(as<double>(regParam["eta"])) ;
    get_lambda_seq(penalty_.lambda_max(data_.XTy_, lambda_factor_), regParam) ;
    
    // Set up the optimizer
    solver_ = SparseOptimizer<matrix,norm>(penalty_, control) ;

  }

template <typename matrix, SparseNorm norm>
double SparseRegularizer<matrix,norm>::get_df() {
  
  double df = set_.size() + data_.centered_ ;
  if (gamma_ > 0) {
    // loop due to sparse encoding. should iterate over the n_zeros only...
    mat SAA(set_.size(),set_.size()) ;
    for (uword i=0;i<set_.size();i++){
      for (uword j=i;j<set_.size();j++){
        SAA(i,j) = data_.S_.at(set_.A_(i),set_.A_(j));
        SAA(j,i) = SAA(i,j);
      }
    }
    // note that XATXAinv_ is in fact (XATXA + lambda S)^-1
    df -= trace(SAA * set_.XATXAinv_); 
  }
  
  return(df);
}

template <typename matrix, SparseNorm norm>
List SparseRegularizer<matrix,norm>::solution_path(const List& control) {
  
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
      nzeros_ = join_cols(nzeros_, beta_/(data_.norm_X_(set_.A_)));
      intercept_.push_back(data_.y_bar_ - dot(beta_, data_.X_bar_(set_.A_))); // X_bar is scaled
      // compute and store debiased coefficients
      set_.inverse_Gram() ;
      beta_debiased_ = set_.XATXAinv_ * (data_.XTy_(set_.A_) - data_.X_bar_(set_.A_) * accu(data_.y_)) ;
      debiased_ = join_cols(debiased_, beta_debiased_/(data_.norm_X_(set_.A_)));
      intercept_debiased_.push_back(data_.y_bar_ - dot(beta_debiased_, data_.X_bar_(set_.A_))) ;
      // store degrees fo freedom and current active set
      df_.push_back(get_df()) ;
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

