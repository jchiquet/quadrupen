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
  
  uword iter = 0, iter_in= 0 ; // count the number of systems solved
  
  vec grad = - data.XTy_(set.A_) + set.XATXA_ * beta0 ;
  vec theta = -sign(grad)   ; // vector of sign of the solution
  
  uvec A = find(abs(beta0) > ZERO) ; // vector of locally active variables
  theta.elem(A)   = sign(beta0.elem(A));
  
  // Solving the quadratic problem
  vec betak = beta0;
  if (set.use_chol_) {
    // betak = solve(trimatu(set.R_), 
    //               solve(trimatl(strans(set.R_)), data.XTy_(set.A_) - lambda * theta));
    betak =  set.Rinv_* set.Rinv_.t() * (data.XTy_(set.A_) - lambda * theta);
  } else {
    iter_in = this->conjugate_gradient(betak, set.XATXA_, data.XTy_(set.A_) - lambda*theta,
                       accuracy, max_iter) ;
  }
  
  // Check for swapping variables
  uvec swap = find(abs(sign(betak.elem(A)) - theta.elem(A)) > ZERO);
  uvec null ;
  if (swap.is_empty()) {
    null = swap ; // this is empty
    beta0 = betak;
  } else {
    swap = A.elem(swap);
    colvec beta0_swap = beta0.elem(swap);
    colvec betak_swap = betak.elem(swap);
    // first, go to zero for the swapped variable which cost the minimum
    vec eta = -beta0_swap / (betak_swap-beta0_swap);
    uword i_swap = eta.index_min();
    double scale = eta(i_swap) ;
    null = swap[i_swap];
    betak = beta0 + (betak-beta0) * scale ;
    // second, solve the problem after swaping the signs of the
    // incriminated variable
    beta0 = betak;
    beta0(null[0]) = -betak_swap[i_swap];
    
    A = find(abs(beta0) > ZERO) ; // new vector of active variables
    theta = -sign(grad)        ; // vector of sign of the solution
    theta.elem(A) = sign(beta0.elem(A));
    
    vec betal = betak;
    if (set.use_chol_) {
      // betal = solve(trimatu(set.R_),
      //               solve( trimatl(strans(set.R_)), data.XTy_(set.A_) - lambda * theta));
      betal =  set.Rinv_* set.Rinv_.t() * (data.XTy_(set.A_) - lambda * theta);
    } else {
      iter_in = this->conjugate_gradient(betal, set.XATXA_, data.XTy_(set.A_) - lambda*theta,
                         accuracy, max_iter) ;
    }
    iter++;
    
    // This is the gradient on the active part of the parameters
    grad = - data.XTy_(set.A_) + set.XATXA_ * betal ;
    // if the sign is coherent, keep that one...
    if (fabs(grad(null[0]) + lambda * as_scalar(sign(betal(null)))) <= ZERO) {
      null = swap; 
      beta0 = betal ;
    } else {
      // otherwise, backtrack to betak
      beta0 = betak ;
    }
  }

  if (set.use_chol_) {
    return(iter) ;
  } else return(iter_in) ;

}

#endif

