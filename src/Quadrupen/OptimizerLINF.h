/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_OPTIMIZER_L1_H
#define _quadrupen_OPTIMIZER_L1_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "GenericOptimizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class OptimizerLINF :
  public GenericOptimizer<matrix,Norm::LINF> {
  
  using GenericOptimizer<matrix,Norm::LINF>::penalty_ ;
  
public:
  
  OptimizerLINF() {} ;
  OptimizerLINF(Penalty<Norm::LINF>&) ;
  
  uword quadratic_breg(
      vec& beta,
      const double& lambda,
      RegressionData<matrix> &data,
      ActiveSet<matrix>& set,
      mat& XTX,
      const double& accuracy,
      const uword& max_iter) ;
  
};

template <typename matrix>
inline OptimizerLINF<matrix>::OptimizerLINF(
    Penalty<Norm::LINF>& penalty) : GenericOptimizer<matrix, Norm::LINF>(penalty) 
    {}

template <typename matrix>
uword OptimizerLINF<matrix>::quadratic_breg(
    vec& beta,
    const double& lambda,
    RegressionData<matrix> &data,
    ActiveSet<matrix>& set,
    mat& XTX,
    const double& accuracy,
    const uword& max_iter) {
  
  uvec B = regspace<uvec>(0,beta.n_elem-1) ; B.shed_rows(set.A_) ;
  vec grad = -data.XTy_ + XTX * beta ;
  vec theta = -sign(grad.elem(B)) ; // sign of the guys on the boundary
  
  uword iter = 0, iter_in= 0 ; // count the number of systems solved
  bool success = false ;
  while ((!success) && (iter < 3)) {
    
    iter++;
    
    // SOLVE THE QUADRATIC PROBLEM
    vec XX_B = XTX.cols(B) * theta;
    if (set.A_.is_empty()) {
      double b  = (dot(theta, data.XTy_.elem(B)) - lambda);
      beta.elem(B) = theta * (b/sum(theta % XX_B.elem(B),0)) ;
    } else {
      // Constructing the system (KKT)
      mat XX = join_rows(
        join_cols(sum(theta % XX_B.elem(B),0), XX_B.elem(set.A_)),
        join_cols(strans(XX_B.elem(set.A_)), set.XATXA_));
      vec b = join_cols(theta.t()*data.XTy_.elem(B)-lambda, data.XTy_.elem(set.A_)) ;
      vec tmp ;
      if (set.use_chol_) {
        // Solving via Cholesky factorization...
        mat R = chol(XX) ;
        tmp = solve(trimatu(R), solve(trimatl(strans(R)), b)) ;
      } else {
        double bound = penalty_.pen_norm(beta.elem(B)) ;
        tmp = join_cols(ones(1) * bound, beta.elem(set.A_));
        iter_in = this->conjugate_gradient(tmp, XX, b, accuracy, max_iter) ;
      }
      beta.elem(B) = theta * tmp[0] ;
      beta.elem(set.A_) = tmp.subvec(1,tmp.n_elem-1) ;
    }
    
    // Handling guys reaching the boundary (leaving the active set)
    double bound = penalty_.pen_norm(beta.elem(B)) ; // current boundary
    uvec ind_toB = find(abs(beta.elem(set.A_)) > bound) ;
    if (!ind_toB.is_empty()) {
      uvec toB = set.A_(ind_toB) ;
      set.del_vars(ind_toB) ;
      B = join_cols(B,toB);
      beta.elem(B) = bound * sign(beta.elem(B));
      theta = sign(beta.elem(B));
    } else {
      success = true ;
    }
  }
  
  // Guys leaving the boundary after optimization (activation)
  grad = -data.XTy_ + XTX * beta ;
  uvec ind_toA  = find(theta == sign(grad.elem(B)));
  if (!ind_toA.is_empty()) {
    uvec toA = B(ind_toA) ;
    set.add_vars(toA, data) ;
    B.shed_rows(ind_toA) ;
    if (B.is_empty()) {
      throw std::runtime_error("Every variable left the boudary. Try more regularization.");
    }
  }
  if (set.use_chol_) {
    return(iter) ;
  } else return(iter_in) ;
  
}

#endif

