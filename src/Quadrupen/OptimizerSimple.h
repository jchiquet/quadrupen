/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_SIMPLE_OPTIMIZER_H
#define _quadrupen_SIMPLE_OPTIMIZER_H

#include "Optimizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, SimpleNorm norm> class SimpleOptimizer: public Optimizer<matrix> {
  
public:
  
  SimpleOptimizer() {} ;
  SimpleOptimizer(SimplePenalty<norm>&, const List&) ;
  
  SimplePenalty<norm> penalty_ ;
  using Optimizer<matrix>::algorithm_  ;
  using Optimizer<matrix>::accuracy_   ;
  using Optimizer<matrix>::maxiter_    ;
  using Optimizer<matrix>::maxfeat_    ;
  using Optimizer<matrix>::verbosity_  ;
  using Optimizer<matrix>::monitoring_ ;
  using Optimizer<matrix>::iter_       ;
  using Optimizer<matrix>::inner_iter_ ;
  using Optimizer<matrix>::gap_        ;
  using Optimizer<matrix>::J_          ;
  using Optimizer<matrix>::D_          ;
  using Optimizer<matrix>::J_vec_      ;
  using Optimizer<matrix>::D_vec_      ;

  using Optimizer<matrix>::optimality_gap ;
  using Optimizer<matrix>::fista ;
  using Optimizer<matrix>::fista_LM ;
  
  uword solve(
      vec& beta,
      vec& grad,
      const double& lambda, 
      const vec& weights,
      const double& gamma, 
      RegressionData<matrix> &data,
      ActiveSet<matrix>& set
  ) ;

  virtual uword quadratic(
      vec &beta,
      vec &grad,
      const double &lambda ,
      RegressionData<matrix> &data,
      ActiveSet<matrix> &set,
      const double& accuracy,
      const uword& max_iter) = 0 ;
  
};

template <typename matrix, SimpleNorm norm>
SimpleOptimizer<matrix, norm>::SimpleOptimizer(
    SimplePenalty<norm>& penalty, const List& control) : 
  Optimizer<matrix>(control) {
    penalty_ = penalty ;
  }

template <typename matrix, SimpleNorm norm>
uword SimpleOptimizer<matrix,norm>::solve(
    vec& beta,
    vec& grad,
    const double& lambda,
    const vec& weights,
    const double& gamma, 
    RegressionData<matrix> &data,
    ActiveSet<matrix>& set) {
  
  if (verbosity_) Rprintf("\n current penalty = %f",lambda) ;
  if (verbosity_) Rprintf("\n nb active variables = %i\n", set.size()) ;

  vec lambda_w = lambda * weights ;
  vec optimality = penalty_.elt_dual_norm(grad) - lambda_w;
  uword var_in = optimality.index_max() ; // highest violation of KKT conditions 
  uword status = 0 ; iter_ = 0 ; bool success = true ; 
  gap_ = std::max(0.0, optimality(var_in)) ;
  J_ = datum::inf ; D_ = datum::inf ;

  auto prox = [this](vec x, vec L) {
    return(penalty_.proximal(x, L));
  } ;
  
  while ((gap_ > accuracy_) && (iter_ <= maxiter_)) {
    R_CheckUserInterrupt();
    iter_++;
    
    // VARIABLE ACTIVATION IF APPLICABLE
    if (set.is_in_[var_in] == 0) { // Is var_in already in the active set?
      set.add_var(var_in, data) ;
      beta.resize(beta.size()+1) ; // update the vector of active parameters
      beta.tail(1) = - 1e-3 * sign(grad(var_in)) ;
      if (verbosity_) {Rprintf("\tnewly added variable %i\n",var_in);}
    }
    
    // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    if (algorithm_ == FISTA) {
      inner_iter_.push_back(
        // fista(beta, grad, lambda, data, set, prox, 1e-10, 10000)
        fista_LM(beta, grad, lambda_w(set.A_), data, set, prox, 1e-10, 10000)
      );
    } else { // QUADRA solver
      try {
        inner_iter_.push_back(
          quadratic(beta, grad, lambda, data, set, 1e-5, 10000)
        );
      } catch (std::runtime_error& error) {
        if (verbosity_ > 0) {
          Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
        }
        success = false ;
      }
    }
    
    // OPTIMALITY TESTING
    grad = - data.XTy_ + set.XTXA_ * beta ;
    optimality = penalty_.elt_dual_norm(grad) - lambda ;
    var_in = optimality.index_max() ;
    gap_ = std::max(0.0, optimality(var_in)) ;
    
    if (monitoring_ > 0) {
      optimality_gap(beta, grad, lambda, gamma, data, set, monitoring_) ;
      J_vec_.push_back(J_) ;
      D_vec_.push_back(D_) ;
    }
    
  }
  if (verbosity_) Rprintf("\tcurrent gap = %f\n",gap_) ;
  
  // Checking convergence status
  if (iter_ >= maxiter_)     { status = 1 ; }
  if (set.size() > maxfeat_) { status = 2 ; }
  if (!success)              { status = 3 ; }
  
  return status ;
}

#endif

