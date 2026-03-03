/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _ElasticNet_H
#define _ElasticNet_H

#include "RegularizerSparse.h"
#include "ActiveSet.h"
#include "OptimizerL1.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class ElasticNet : 
  public SimpleSparseRegularizer<matrix,SimpleNorm::L1>{
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
    using SparseRegularizer<matrix>::intercept_debiased_   ;
    using SparseRegularizer<matrix>::iA_   ;
    using SparseRegularizer<matrix>::jA_   ;
    using SimpleSparseRegularizer<matrix,SimpleNorm::L1>::penalty_   ;
    using SimpleSparseRegularizer<matrix,SimpleNorm::L1>::set_ ;
    using SimpleSparseRegularizer<matrix,SimpleNorm::L1>::get_lambda_max ;
    
  ElasticNet(const RegressionData<matrix>&, const List&, const List&);

  // Specific to Elastic-Net regularization
  OptimizerL1<matrix> solver_ ; // Solvers for L1 penalty

  List solution_path(const List&);

  void optimality_gap(double lambda_, uword type) ;
  
  // Compute degrees of freedom for the current estimate
  double get_df() ;

};

template <typename matrix>
ElasticNet<matrix>::ElasticNet(
  const RegressionData<matrix>& data, const List& regParam, const List& control) :
  SimpleSparseRegularizer<matrix,SimpleNorm::L1>::SimpleSparseRegularizer(data, regParam) {
    
    // set the penalty to l1
    penalty_ = SimplePenalty<SimpleNorm::L1>() ;
    get_lambda_seq(get_lambda_max(), regParam) ;

    // Set up the optimizer
    solver_ = OptimizerL1<matrix>(penalty_, control) ;
    
    // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
    data_.scale_struct(sqrt(gamma_)*pow(lambda_factor_,-1)) ;
    data_.scale_regressors(lambda_factor_) ;
    
    // Initialize the active set, beta_ and gradient with starting coefficient
    vec beta0 = control["beta0"] ;
    uvec A0 = find(beta0) ;
    if (A0.is_empty()) {
      set_  = ActiveSet(data_, as<bool>(control["usechol"])) ;
      grad_ = - data_.XTy_ ;
    } else {
      set_  = ActiveSet(data_, A0, as<bool>(control["usechol"])) ;
      beta_ = beta0(A0) ;
      grad_ = - data_.XTy_ + set_.XTXA_ * beta_  ;
    }
  }

template <typename matrix>
double ElasticNet<matrix>::get_df() {
  
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

template <typename matrix>
List ElasticNet<matrix>::solution_path(const List& control) {
  
  // Variables monitoring the algorithm
  vector<double> gap, timing ; // timings and optimality measures
  vector<uword> status, iactive ; // convergence and # of inner/outer iterates

  // LAMBDA LOOP
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda_ : lambdas_) {
    if (solver_.verbosity_) {
      Rprintf("\n lambda_l1 = %f",lambda_) ;
      Rprintf("\n nb active variables = %i\n", set_.size()) ;
    }
    
    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    status.push_back(
      solver_.solve(beta_, grad_, lambda_, gamma_, data_, set_)
    ) ;
    gap.push_back(solver_.gap_) ;
    iactive.push_back(solver_.iter_) ;
    
    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      set_.inverse_Gram() ;
      nzeros_ = join_cols(nzeros_, beta_/(data_.norm_X_(set_.A_) % lambda_factor_(set_.A_)));
      vec beta_debiased = set_.XATXAinv_ * (data_.XTy_(set_.A_) - data_.X_bar_(set_.A_) * accu(data_.y_)) ;
      debiased_ = join_cols(debiased_, beta_debiased/(data_.norm_X_(set_.A_) % lambda_factor_(set_.A_)));
      intercept_.push_back(data_.y_bar_ - dot(beta_, data_.X_bar_(set_.A_)));
      intercept_debiased_.push_back(data_.y_bar_ - dot(beta_debiased, data_.X_bar_(set_.A_))) ;
      iA_ = join_rows(iA_, set_.A_.t()) ;
      jA_ = join_rows(jA_, df_.size()*ones<urowvec>(set_.size()) );
      df_.push_back(get_df()) ;
    }
    
    timing.push_back(timer.toc()) ;
  } // END OF THE LOOP OVER LAMBDA
  
  lambdas_.resize(df_.size()) ;
  
  return(
    List::create(
      Named("it_active")      = iactive,
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

