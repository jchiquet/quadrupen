/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_GROUP_OPTIMIZER_H
#define _quadrupen_GROUP_OPTIMIZER_H

#include "Optimizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, MixedNorm norm> 
class GroupOptimizer : public Optimizer<matrix >{
public:
  
  GroupOptimizer() {} ;
  GroupOptimizer(MixedPenalty<norm>&) ;
  
  uword fista(
      vec& beta,
      const double& lambda, 
      RegressionData<matrix> &data,
      ActiveSetGroup<matrix>& set,
      const double& accuracy, 
      const uword& max_iter) ;
  
protected:
  MixedPenalty<norm> penalty_ ;
  
};

template <typename matrix, MixedNorm norm>
inline GroupOptimizer<matrix, norm>::GroupOptimizer(
    MixedPenalty<norm>& penalty) : 
  penalty_ (penalty), Optimizer<matrix>() {}

template <typename matrix, MixedNorm norm>
uword GroupOptimizer<matrix,norm>::fista(
    vec& beta,
    const double& lambda,
    RegressionData<matrix> &data,
    ActiveSetGroup<matrix>& set,
    const double& accuracy,
    const uword& max_iter) {
  
  vec betak = beta  ; // output vector
  vec betal = beta  ;
  double delta = 2*accuracy  ; // change in beta
  double L = max( set.XATXA_.diag()) ; // Lipchitz constant
  
  double t0 = 1.0, tk ; // auxiliary variables in FISTA 
  uword iter = 0      ; // current iterate
  while ((delta > accuracy/beta.n_elem ) && (iter < max_iter)) {
    
    double l_num, l_den ;
    double f0, fk ;
    vec XATXA_betal = set.XATXA_ * betal ;
    f0 = dot(betal, .5 * XATXA_betal  - data.XTy_(set.A_)) ;
    vec grad = XATXA_betal - data.XTy_(set.A_);
    
    // Line search over L
    bool found=false;
    while(!found) {
      // Apply proximal operator (implemented in penalty object)
      vec prox_arg = betal - grad/L;
      betak = penalty_.proximal(prox_arg, lambda/L, set.grp_sizes_(set.G_));
      
      fk = dot(betak, .5 * set.XATXA_ * betak - data.XTy_(set.A_)) ;
      l_num = 2 * (fk - f0 - dot(grad, betak-betal));
      l_den = accu(pow(betak-betal,2));
      
      if ((L * l_den >= l_num) || (sqrt(l_den) < accuracy)) {
        found = true;
      } else {
        L = fmax(2*L, l_num/l_den);
      }
      
      R_CheckUserInterrupt();
    }
    
    // updating t
    tk = 0.5 * (1+sqrt(1+4*t0*t0));
    
    // updating s
    betal = betak + (t0-1)/tk * ( betak - beta );
    
    // preparing next iterate
    delta = sqrt(l_num);
    beta = betak;
    t0 = tk;
    iter++;
    
    R_CheckUserInterrupt();
  }
  
  return(iter) ;
  
}

#endif

