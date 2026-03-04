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

template <typename matrix, MixedNorm norm> class GroupOptimizer : public Optimizer<matrix >{
  
public:
  
  GroupOptimizer() {} ;
  GroupOptimizer(MixedPenalty<norm>&, const List&) ;
  
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
  MixedPenalty<norm> penalty_ ;

  using Optimizer<matrix>::optimality_gap ;
  
  uword solve(
      vec& beta,
      vec& grad,
      const double& lambda,
      const double& gamma, 
      RegressionData<matrix> &data,
      ActiveSetGroup<matrix>& set
  ) ;
  
  uword fista(
      vec& beta,
      vec& grad,
      const double& lambda, 
      RegressionData<matrix> &data,
      ActiveSetGroup<matrix>& set,
      const double& accuracy, 
      const uword& max_iter
  ) ;

};

template <typename matrix, MixedNorm norm>
GroupOptimizer<matrix, norm>::GroupOptimizer(
    MixedPenalty<norm>& penalty, const List& control) : 
  penalty_ (penalty), Optimizer<matrix>(control) {}

template <typename matrix, MixedNorm norm>
uword GroupOptimizer<matrix,norm>::solve(
    vec& beta,
    vec& grad,
    const double& lambda,
    const double& gamma, 
    RegressionData<matrix> &data,
    ActiveSetGroup<matrix>& set) {
  
  if (verbosity_) Rprintf("\n current penalty = %f",lambda) ;
  if (verbosity_) Rprintf("\n nb active groups = %i\n", set.size_grp()) ;
  
  vec optimality = penalty_.elt_norm(grad, set.grp_sizes_) - lambda;
  uword grp_in = optimality.index_max() ; // highest violation of KKT conditions 
  uword status = 0 ; iter_ = 0 ; bool success = true ; 
  gap_ = std::max(0.0, optimality(grp_in)) ;
  J_ = datum::inf ; D_ = datum::inf ;
  
  while ((gap_ > accuracy_) && (iter_ <= maxiter_)) {
    R_CheckUserInterrupt();
    iter_++;
    
    // VARIABLE ACTIVATION IF APPLICABLE
    if (set.is_grp_in_[grp_in] == 0) { // Is var_in already in the active set?
      set.add_group(grp_in, data) ;
      beta.resize(beta.size()+set.grp_sizes_(grp_in)) ; // update the vector of active parameters
      beta.tail(set.grp_sizes_(grp_in)) = - 1e-3 * sign(grad(set.group_[grp_in])) ;
      if (verbosity_) {Rprintf("\tnewly added group %i\n",grp_in);}
    } else if (verbosity_) {Rprintf("\talready in %i\n",grp_in);}
    
    // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    // if (algorithm_ == FISTA) {
      inner_iter_.push_back(
        fista(beta, grad, lambda, data, set, 1e-10, 10000)
      );
    // }
    
    // OPTIMALITY TESTING
    grad = - data.XTy_ + set.XTXA_ * beta ;
    optimality = penalty_.elt_norm(grad, set.grp_sizes_) - lambda ;
    grp_in = optimality.index_max() ;
    gap_ = std::max(0.0, optimality(grp_in)) ;
    
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

template <typename matrix, MixedNorm norm>
uword GroupOptimizer<matrix,norm>::fista(
    vec& beta,
    vec& grad,
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
    grad(set.A_) = - data.XTy_(set.A_) + XATXA_betal ;
    
    // Line search over L
    bool found=false;
    while(!found) {
      // Apply proximal operator (implemented in penalty object)
      vec prox_arg = betal - grad(set.A_)/L;
      betak = penalty_.proximal(prox_arg, lambda/L, set.grp_sizes_(set.G_));
      
      fk = dot(betak, .5 * set.XATXA_ * betak - data.XTy_(set.A_)) ;
      l_num = 2 * (fk - f0 - dot(grad(set.A_), betak-betal));
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

