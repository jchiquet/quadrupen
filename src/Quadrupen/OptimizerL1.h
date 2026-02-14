/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_OPTIMIZER_L1_H
#define _quadrupen_OPTIMIZER_L1_H

#include "GenericOptimizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class OptimizerL1 :
  public GenericOptimizer<matrix,Norm::L1> {

  using GenericOptimizer<matrix,Norm::L1>::penalty_ ;
  
  public:
  
  OptimizerL1() {} ;
  OptimizerL1(Penalty<Norm::L1>&) ;

  uword quadratic_enet(
      vec &beta0,
      const double &lambda ,
      RegressionData<matrix> &data,
      ActiveSet<matrix> &set,
      const double& accuracy,
      const uword& max_iter) ;

};

template <typename matrix>
inline OptimizerL1<matrix>::OptimizerL1(
    Penalty<Norm::L1>& penalty) : GenericOptimizer<matrix, Norm::L1>(penalty) 
  {}

template <typename matrix>
uword OptimizerL1<matrix>::quadratic_enet(
    vec &beta0,
    const double &lambda ,
    RegressionData<matrix> &data,
    ActiveSet<matrix> &set,
    const double& accuracy,
    const uword& max_iter) {
  
  uword iter = 1, iter_in= 0 ; // count the number of systems solved
  
  // Solving the quadratic problem
  vec betak = beta0;
  vec theta = sign(beta0) ; // vector of sign of the solution
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
    beta0 = betak;
  } else {
    iter++;
    vec beta0_swap = beta0(swap); vec betak_swap = betak(swap);
    
    // first, go to zero for the swapped variable which cost the minimum
    vec eta = -beta0_swap / (betak_swap-beta0_swap);
    uword i_min = eta.index_min();
    betak = beta0 + (betak-beta0) * eta(i_min) ;
    
    // second, solve the problem after swapping the signs of the incriminated variable
    uword i_swap = swap[i_min];
    beta0 = betak;
    beta0(i_swap) = -betak_swap[i_min];
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
        beta0 = betal ;
    } else {
      beta0 = betak ; // otherwise, backtrack to betak
      // if (verbose) Rprint("\tremoving variables %i", set.A_(i_swap)) ;
      set.del_var(i_swap) ; // and desactivate de zeroed variable
      beta0.shed_row(i_swap) ;
    }
  }

  if (set.use_chol_) {
    return(iter) ;
  } else return(iter_in) ;

}

#endif

