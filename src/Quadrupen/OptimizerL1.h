/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_OPTIMIZER_L1_H
#define _quadrupen_OPTIMIZER_L1_H

#include "Optimizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class OptimizerL1 :
  public SimpleOptimizer<matrix,SimpleNorm::L1> {

  using SimpleOptimizer<matrix,SimpleNorm::L1>::penalty_ ;
  
  public:
  
  OptimizerL1() {} ;
  OptimizerL1(SimplePenalty<SimpleNorm::L1>&) ;
  OptimizerL1(SimplePenalty<SimpleNorm::L1>&, const List& control) ;
  
  uword quadratic(
      vec &beta,
      vec &grad,
      const double &lambda ,
      RegressionData<matrix> &data,
      ActiveSet<matrix> &set,
      const double& accuracy,
      const uword& max_iter) override ;

};

template <typename matrix>
OptimizerL1<matrix>::OptimizerL1(
    SimplePenalty<SimpleNorm::L1>& penalty) : SimpleOptimizer<matrix, SimpleNorm::L1>(penalty) 
  {}

template <typename matrix>
OptimizerL1<matrix>::OptimizerL1(SimplePenalty<SimpleNorm::L1>& penalty, const List& control) : 
  SimpleOptimizer<matrix, SimpleNorm::L1>(penalty, control) 
  {}

template <typename matrix>
uword OptimizerL1<matrix>::quadratic(
    vec &beta,
    vec &grad,
    const double &lambda ,
    RegressionData<matrix> &data,
    ActiveSet<matrix> &set,
    const double& accuracy,
    const uword& max_iter) {
  
  uword iter = 1, iter_in= 0 ; // count the number of systems solved
  
  // Solving the quadratic problem
  vec betak = beta;
  vec theta = sign(beta) ; // vector of sign of the solution
  if (set.use_chol_) {
    betak =  set.Rinv_* set.Rinv_.t() * (data.XTy_(set.A_) - lambda * theta);
  } else {
    iter_in = 
      this->conjugate_gradient(betak, set.XATXA_, 
                               data.XTy_(set.A_) - lambda*theta,
                               accuracy, max_iter) ;
  }
  
  // Check for swapping variables
  uvec swap = find(abs(sign(betak) - theta) > ZERO);
  if (swap.is_empty()) {
    beta = betak;
  } else {
    iter++;
    vec beta_swap = beta(swap); vec betak_swap = betak(swap);
    
    // first, go to zero for the swapped variable which cost the minimum
    vec eta = -beta_swap / (betak_swap-beta_swap);
    uword i_min = eta.index_min();
    betak = beta + (betak-beta) * eta(i_min) ;
    
    // second, solve the problem after swapping the signs of the incriminated variable
    uword i_swap = swap[i_min];
    beta = betak;
    beta(i_swap) = -betak_swap[i_min];
    double grad_swap = - data.XTy_(set.A_(i_swap)) + 
      dot(set.XATXA_.row(i_swap), betak) ;
    theta(i_swap) = - sign(grad_swap) ;
    vec betal = betak;
    if (set.use_chol_) {
      betal =  set.Rinv_* set.Rinv_.t() * (data.XTy_(set.A_) - lambda * theta);
    } else {
      iter_in = 
        this->conjugate_gradient(betal, set.XATXA_, 
                                 data.XTy_(set.A_) - lambda*theta,
                                 accuracy, max_iter) ;
    }

    // if the sign is coherent, keep that one...
    grad_swap = - data.XTy_(set.A_(i_swap)) + 
      dot(set.XATXA_.row(i_swap), betal) ;
    if (abs(grad_swap + lambda * sign(betal(i_swap))) <= ZERO) {
        beta = betal ;
    } else {
      beta = betak ; // otherwise, backtrack to betak
      // if (verbose) Rprint("\tremoving variables %i", set.A_(i_swap)) ;
      set.del_var(i_swap) ; // and desactivate de zeroed variable
      beta.shed_row(i_swap) ;
    }
  }
  
  if (set.use_chol_) {
    return(iter) ;
  } else return(iter_in) ;

}

#endif

