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
  using Optimizer<matrix>::fista ;
  using Optimizer<matrix>::fista_LM ;
  
  uword solve(
      vec& beta,
      vec& grad,
      const double& lambda,
      const vec& weights,
      const double& gamma, 
      RegressionData<matrix> &data,
      ActiveSetGroup<matrix>& set
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
    const vec& weights,
    const double& gamma, 
    RegressionData<matrix> &data,
    ActiveSetGroup<matrix>& set) {
  
  if (verbosity_) Rprintf("\n current penalty = %f",lambda) ;
  if (verbosity_) Rprintf("\n nb active groups = %i\n", set.size_grp()) ;
  
  vec lambda_w = lambda * weights ;
  vec optimality = penalty_.elt_dual_norm(grad, set.grp_sizes_) - lambda_w;
  uword grp_in = optimality.index_max() ; // highest violation of KKT conditions 
  uword status = 0 ; iter_ = 0 ; bool success = true ; 
  gap_ = std::max(0.0, optimality(grp_in)) ;
  J_ = datum::inf ; D_ = datum::inf ;
  
  while ((gap_ > accuracy_) && (iter_ <= maxiter_)) {
    R_CheckUserInterrupt();
    iter_++;
    
    // VARIABLE ACTIVATION IF APPLICABLE
    if (set.is_grp_in_[grp_in] == 0 && optimality(grp_in) > 0) { // Is var_in already in the active set?
      set.add_group(grp_in, data) ;
      beta.resize(beta.size()+set.grp_sizes_(grp_in)) ; // update the vector of active parameters
      beta.tail(set.grp_sizes_(grp_in)) = - 1e-3 * sign(grad(set.group_[grp_in])) ;
      if (verbosity_) {Rprintf("\tnewly added group %i\n",grp_in);}
    } 
    
    // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    // if (algorithm_ == FISTA) {
    auto prox = [this, set](vec x, vec L) {
      return(penalty_.proximal(x, L, set.grp_sizes_(set.G_)));
    } ;
    inner_iter_.push_back(
      fista_LM(beta, grad, lambda_w(set.G_), data, set, prox, 1e-6, 10000)
    );
    // }

    // VARIABLE DELETION IF APPLICABLE
    uvec vanish = find(
      penalty_.elt_norm(grad(set.A_), set.grp_sizes_(set.G_)) < lambda_w(set.G_) + ZERO &&
      penalty_.elt_norm(beta, set.grp_sizes_(set.G_)) < ZERO
    ) ;
    if (!vanish.is_empty()) { // Is var_in already in the active set?
      set.del_group(vanish(0), beta) ;
      if (verbosity_) {Rprintf("\tremoved group %i\n",set.G_(vanish(0)));}
    } 

    // OPTIMALITY TESTING
    grad = - data.XTy_ + set.XTXA_ * beta ;
    optimality = penalty_.elt_dual_norm(grad, set.grp_sizes_) - lambda_w ;
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

#endif

