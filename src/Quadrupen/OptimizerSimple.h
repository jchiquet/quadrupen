/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_SIMPLE_OPTIMIZER_H
#define _quadrupen_SIMPLE_OPTIMIZER_H

#include "Optimizer.h"
#include "PenaltySimple.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, SimpleNorm norm> class SimpleOptimizer: public Optimizer {
  
public:
  
  SimpleOptimizer() {} ;
  SimpleOptimizer(SimplePenalty<norm>&, const List&) ;
  
  SimplePenalty<norm> penalty_ ;
  using Optimizer::algorithm_  ;
  using Optimizer::accuracy_   ;
  using Optimizer::maxiter_    ;
  using Optimizer::maxfeat_    ;
  using Optimizer::verbosity_  ;
  using Optimizer::monitoring_ ;
  using Optimizer::iter_       ;
  using Optimizer::inner_iter_ ;
  using Optimizer::gap_        ;
  using Optimizer::J_          ;
  using Optimizer::D_          ;
  using Optimizer::J_vec_      ;
  using Optimizer::D_vec_      ;
  using Optimizer::optimality_gap ;
  using Optimizer::fista ;
  using Optimizer::pgd ;
  
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
      const vec &lambda,
      const vec &XTy,
      ActiveSet<matrix> &set,
      const double& accuracy,
      const uword& max_iter) = 0 ;
  
};

template <typename matrix, SimpleNorm norm>
SimpleOptimizer<matrix, norm>::SimpleOptimizer(
    SimplePenalty<norm>& penalty, const List& control) : 
    Optimizer(control) {
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

  vec optimality = penalty_.optimality(grad, lambda, weights) ;
  uword var_in = optimality.index_max() ; // highest violation of KKT conditions 
  uword status = 0 ; iter_ = 0 ; bool success = true ; 
  gap_ = std::max(0.0, optimality(var_in)) ;
  J_ = datum::inf ; D_ = datum::inf ;

  while ((gap_ > accuracy_) && (iter_ <= maxiter_)) {
    R_CheckUserInterrupt();
    iter_++;

    // VARIABLE ACTIVATION IF APPLICABLE
    if (set.is_in_[var_in] == 0) { // Is var_in already in the active set?
      set.add_var(var_in, data) ;
      beta.insert_rows(beta.n_elem, 1); // update the vector of active parameters
      if (algorithm_ ==  QUADRA) {
        beta.tail(1).fill(- 1e-3 * sign(grad(var_in)));
      } else {
        beta.tail(1).fill(0.0);
      }
      if (verbosity_) {Rprintf("\tnewly added variable %i\n",var_in);}
    }
    
    // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    if (algorithm_ ==  QUADRA) { // Newton-based solver
      inner_iter_.push_back(
        quadratic(beta, lambda*weights, data.XTy_, set, 1e-9, 1000)
      );
      grad = - data.XTy_ + set.XTXA_ * beta ;
    } 
    else { // Proximal-based solvers
      auto prox = [this, &set, &weights](const vec& x, double l) {
        return(penalty_.proximal(x, l, weights.elem(set.A_)));
      };
      vec beta_old = beta ;
      if (algorithm_ == FISTA) {
        inner_iter_.push_back(
          fista(beta, lambda, data.XTy_.elem(set.A_), set.XATXA_, prox, 1e-7, 3000)
        );
      } else if (algorithm_ == PGD) {
        inner_iter_.push_back(
          pgd(beta, lambda, data.XTy_.elem(set.A_), set.XATXA_, prox, 1e-7, 3000, 5)
        );
      }
      grad += set.XTXA_ * (beta - beta_old); // Incremental update of the gradient
      uvec vanish = find( // Variable deletion if applicable
        penalty_.optimality(grad.elem(set.A_), lambda, weights.elem(set.A_)) <= accuracy_ &&
          abs(beta) < accuracy_/10 * weights.elem(set.A_)
      ) ;
      if (!vanish.is_empty()) {
        if (verbosity_) {set.A_(vanish).t().print("Removing variables");}
        set.del_vars(vanish, beta) ;
      }
    }

    // OPTIMALITY TESTING
    optimality = penalty_.optimality(grad, lambda, weights) ;
    var_in = optimality.index_max() ;
    gap_ = std::max(0.0, optimality(var_in)) ;
    
    if (monitoring_ > 0) {
      optimality_gap(beta, grad, lambda, gamma, data.XTy_(set.A_), set.XATXA_, data.norm_y_, set.A_, monitoring_) ;
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

