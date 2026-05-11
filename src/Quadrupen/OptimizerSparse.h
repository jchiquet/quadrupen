/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "Optimizer.h"
#include "PenaltySparse.h"
#include "ActiveSet.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, SparseNorm norm>
class SparseOptimizer: public Optimizer {
  
public:
  
  SparseOptimizer() {} ;
  SparseOptimizer(SparsePenalty<norm>&, const List&) ;
  
  SparsePenalty<norm> penalty_ ;
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
  
  uword working_set(
      vec& beta,
      vec& grad,
      const double& lambda, 
      const vec& weights,
      const double& gamma, 
      RegressionData<matrix> &data,
      ActiveSet<matrix>& set
  ) ;

  uword quadratic(
      vec &beta,
      const vec &lambda,
      const vec &XTy,
      ActiveSet<matrix> &set,
      const double& accuracy,
      const uword& max_iter) ;
  
};

template <typename matrix, SparseNorm norm>
SparseOptimizer<matrix, norm>::SparseOptimizer(
    SparsePenalty<norm>& penalty, const List& control) : 
    Optimizer(control) {
    penalty_ = penalty ;
  }

template <typename matrix, SparseNorm norm>
uword SparseOptimizer<matrix, norm>::quadratic(
    vec &beta,
    const vec &lambda,
    const vec &XTy,
    ActiveSet<matrix> &set,
    const double& accuracy,
    const uword& max_iter) {
  
  uword iter = 0, iter_in= 0 ;
  bool signs_stable = false;
  
  while (!signs_stable && iter < 10 && set.size() > 0) { // Max 10 swaps
    iter++;  
    // Solving the quadratic problem
    vec betak = beta;
    vec theta = sign(beta) ; // vector of sign of the solution
    if (set.use_chol_) {
      // Step 1 - Forward Substitution
      vec tmp = arma::solve(trimatl(set.R_.t()), XTy(set.A_) - lambda(set.A_) % theta);
      // Step 2 - Backward Substitution
      betak = arma::solve(trimatu(set.R_), tmp);
    } else {
      iter_in = 
        this->conjugate_gradient(betak, set.XATXA_, 
                                 XTy(set.A_) - lambda(set.A_) % theta,
                                 accuracy, max_iter) ;
    }
    
    // Check for swapping variables
    // uvec swap = find(abs(betak) > accuracy/100 && sign(betak) != theta) ;
    uvec swap = find(sign(betak) != theta) ;
    if (swap.is_empty()) {
      beta = betak;
      signs_stable = true;
    } else {
      // Find the first variable hitting zero
      vec eta = -beta(swap) / (betak(swap) - beta(swap));
      uword i_min = eta.index_min();
      uword idx_to_remove = swap[i_min];
      
      // Interpolate by moving all variables to this point
      beta = beta + eta(i_min) * (betak - beta);
      
      // Remove the incriminated variable
      if (verbosity_) Rprintf("\tremoving variables %i\n", set.A_(idx_to_remove)) ;
      set.del_var(idx_to_remove, beta) ; 
      
      // Update theta on the new set
      theta = sign(beta);      
    }
  }
  
  if (set.use_chol_) {
    return(iter) ;
  } else return(iter_in) ;
  
}

template <typename matrix, SparseNorm norm>
uword SparseOptimizer<matrix,norm>::working_set(
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


