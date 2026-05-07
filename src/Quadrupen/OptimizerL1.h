/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "OptimizerSparse.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class OptimizerL1 : public SimpleOptimizer<matrix,SparseNorm::L1> {
  
public:
  
  using SimpleOptimizer<matrix,SparseNorm::L1>::penalty_ ;
  using Optimizer::verbosity_ ;
    
  OptimizerL1() {} ;
  OptimizerL1(SparsePenalty<SparseNorm::L1>&, const List& control) ;
  
  uword quadratic(
      vec &beta,
      const vec &lambda,
      const vec &XTy,
      ActiveSet<matrix> &set,
      const double& accuracy,
      const uword& max_iter) override ;
  
};

template <typename matrix>
OptimizerL1<matrix>::OptimizerL1(SparsePenalty<SparseNorm::L1>& penalty, const List& control) : 
  SimpleOptimizer<matrix, SparseNorm::L1>(penalty, control) 
  {}

template <typename matrix>
uword OptimizerL1<matrix>::quadratic(
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
      vec tmp = solve(trimatl(set.R_.t()), XTy(set.A_) - lambda(set.A_) % theta);
      // Step 2 - Backward Substitution
      betak = solve(trimatu(set.R_), tmp);
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


