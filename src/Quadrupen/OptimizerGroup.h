/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_GROUP_OPTIMIZER_H
#define _quadrupen_GROUP_OPTIMIZER_H

#include "Optimizer.h"
#include "PenaltyGroup.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, GroupNorm norm> class GroupOptimizer : public Optimizer {
  
public:
  
  GroupOptimizer() {} ;
  GroupOptimizer(GroupPenalty<norm>&, const List&) ;
  
  GroupPenalty<norm> penalty_ ;
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
      ActiveSetGroup<matrix>& set
  ) ;

  uword quadratic(
      vec &beta,
      const double lambda,
      const vec& weights,
      const vec &XTy,
      ActiveSetGroup<matrix> &set,
      const double& accuracy)  ;

};

template <typename matrix, GroupNorm norm>
GroupOptimizer<matrix, norm>::GroupOptimizer(
    GroupPenalty<norm>& penalty, const List& control) : 
  Optimizer(control) {
  penalty_  = penalty ;
}

template <typename matrix, GroupNorm norm>
uword GroupOptimizer<matrix, norm>::quadratic(
    vec &beta,
    const double lambda,
    const vec& weights,
    const vec &XTy,
    ActiveSetGroup<matrix> &set,
    const double& tol) {
  
  double l1 = lambda * penalty_.alpha_ ;
  double l2 = lambda * (1-penalty_.alpha_) ;
  
  uword iter = 0;
  bool stable = false;
  const uword nb_active_groups = set.size_grp();
  
  while (!stable && iter < 15) {
    iter++;
    vec beta_old = beta;
    uword offset = 0; // go through active variables
    uvec groups_to_remove;
    
    for (uword k = 0; k < nb_active_groups; ++k) {
      uword sz   = set.grp_sizes_(set.G_(k)); // group size
      uvec ind_g = regspace<uvec>(offset, offset + sz - 1);
      
      vec beta_g = beta.subvec(offset, offset + sz - 1);
      double norm_g = std::sqrt(arma::dot(beta_g, beta_g) + 1e-12);
      
      // Only optimize if current group is nonzero
      if (norm_g > 1e-15) {
        
        // Get current group residuals:  XTy_g - X'X_g,A * beta_A
        vec res_g = XTy(ind_g) - set.XATXA_.rows(ind_g) * beta;
        res_g += set.XATXA_(ind_g, ind_g) * beta_g; 
        
        // Newton-Raphson Hessian: H = Xg'Xg + (l2 * wk / ||bg||) * I
        mat Hg = set.XATXA_(ind_g, ind_g);
        Hg.diag() += (l2 * weights(k)) / norm_g;
        
        // Solve Hg * b_new = r_g - l1 * wk * sign(b_g)
        // Add univariate l1 penalty (Sparse Group)
        vec beta_new = arma::solve(Hg, res_g, solve_opts::fast);
        
        // Application du Soft-Thresholding Lasso composante par composante
        if (penalty_.alpha_ > 0) {
          for (uword i = 0; i < sz; ++i) {
            double h_ii = Hg(i,i);
            beta_new(i) = sign(beta_new(i)) * std::max(0.0, std::abs(beta_new(i)) - (l1 * weights(k)) / h_ii);
          }
        }
        
        if (arma::norm(beta_new, 2) < 1e-15) {
          beta_new.zeros();
          groups_to_remove.insert_rows(groups_to_remove.n_elem, 1);
          groups_to_remove.tail(1) = k; // 
        }
        
        beta.subvec(offset, offset + sz - 1) = beta_new;
      }
      
      offset += sz; // go to next group
    }

    if (!groups_to_remove.is_empty()) {
      if (verbosity_) set.G_(groups_to_remove).print("\tremoved group %i\n") ;
      set.del_groups(groups_to_remove, beta);
      // Après suppression, les offsets et la taille de beta ont changé
      // Let working set handle new active set
      return iter; 
    }
    
    if (arma::norm(beta - beta_old, "inf") < tol) stable = true;
    
  }
  return iter;
}

template <typename matrix, GroupNorm norm>
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
  
  vec optimality = penalty_.optimality(grad, lambda, set.grp_sizes_, weights) ;
  uword grp_in = optimality.index_max() ; // highest violation of KKT conditions 
  uword status = 0 ; iter_ = 0 ; bool success = true ; 
  gap_ = std::max(0.0, optimality(grp_in)) ;
  J_ = datum::inf ; D_ = datum::inf ;
  
  while ((gap_ > accuracy_) && (iter_ <= maxiter_)) {
    R_CheckUserInterrupt();
    iter_++;
    double current_tol = 1e-7;

    // VARIABLE ACTIVATION IF APPLICABLE
    if (set.is_grp_in_[grp_in] == 0 && optimality(grp_in) > 0) { // Is var_in already in the active set?
      set.add_group(grp_in, data) ;
      beta.insert_rows(beta.n_elem, set.grp_sizes_(grp_in)); // update the vector of active parameters
      if (algorithm_ ==  QUADRA) {
        beta.tail(set.grp_sizes_(grp_in)).fill(- 1e-3 * sign(grad(grp_in)));
      } else {
        beta.tail(set.grp_sizes_(grp_in)).fill(0.0);
      }
      if (verbosity_) {Rprintf("\tnewly added group %i\n",grp_in);}
    } 

    // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    if (algorithm_ == QUADRA) {
      inner_iter_.push_back(
        quadratic(beta, lambda, weights(set.G_), data.XTy_(set.A_), set, current_tol)
      );
      grad = - data.XTy_ + set.XTXA_ * beta ;
    } else {
      auto prox = [this, &set, &weights](const vec& x, const double l) {
        return(penalty_.proximal(x, l, set.grp_sizes_(set.G_), weights(set.G_)));
      } ;
      vec beta_old = beta ;
      if (algorithm_ == FISTA) {
        inner_iter_.push_back(
          fista(beta, lambda, data.XTy_(set.A_), set.XATXA_, prox, current_tol, 3000)
        );
      } else if (algorithm_ == PGD) {
        inner_iter_.push_back(
          pgd(beta, lambda, data.XTy_(set.A_), set.XATXA_, prox, current_tol, 3000, 3)
        );
      }
      grad += set.XTXA_ * (beta - beta_old);
    }

    // VARIABLE DELETION IF APPLICABLE
    uvec vanish = find(
      penalty_.optimality(grad(set.A_), lambda, set.grp_sizes_(set.G_), weights(set.G_)) <= accuracy_ &&
        penalty_.elt_norm(beta, set.grp_sizes_(set.G_), ones(set.size_grp())) <= accuracy_/10
    ) ;
    if (!vanish.is_empty()) {
      if (verbosity_) set.G_(vanish).print("\tremoved group %i\n") ;
      set.del_groups(vanish, beta) ;
    } 
    
    // OPTIMALITY TESTING
    optimality = penalty_.optimality(grad, lambda, set.grp_sizes_, weights) ;
    grp_in = optimality.index_max() ;
    gap_ = std::max(0.0, optimality(grp_in)) ;
    
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

