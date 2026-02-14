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

  mat updateCholeskyFromExisting(const mat& R, const vec& b) ;
  
};

template <typename matrix>
mat OptimizerLINF<matrix>::updateCholeskyFromExisting(const mat& R, const vec& b) {
  uword p = R.n_cols + 1;
  mat R_ ;
  colvec rp  = zeros<colvec>(p,1);
  rp.subvec(0,p-2) = solve (trimatl(strans(R)), b.subvec(0,p-2));
  rp(p-1) = sqrt(b(p-1) - dot(rp,rp));
  R_ = join_rows( join_cols(R_, zeros<mat>(1,p-1)) , rp);
  return(R_) ;
}

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
  vec theta = -sign(grad(B)) ; // sign of the guys on the boundary
  
  uword iter = 0, iter_in= 0 ; // count the number of systems solved
  bool success = false ;
  while ((!success) && (iter < 3)) {
    
    iter++;
    
    // SOLVE THE QUADRATIC PROBLEM
    vec XX_B = XTX.cols(B) * theta;
    if (set.A_.is_empty()) {
      double b  = (dot(theta, data.XTy_(B)) - lambda);
      beta(B) = theta * (b/sum(theta % XX_B(B),0)) ;
    } else {
      // Constructing the system (KKT)
      vec tmp = join_cols(XX_B(set.A_),sum(theta % XX_B(B),0)) ;
      mat XX = join_rows(join_cols(set.XATXA_, trans(XX_B(set.A_))), tmp);
      vec b = data.XTy_(set.A_) ; b.resize(b.n_elem + 1) ;
      b.tail(1) = dot(theta, data.XTy_(B))-lambda;
      if (set.use_chol_) {
        // Solving via Cholesky factorization...
        mat R = updateCholeskyFromExisting(set.R_, tmp) ;
        tmp = solve(trimatu(R), solve(trimatl(strans(R)), b)) ;
      } else {
        double bound = penalty_.pen_norm(beta(B)) ;
        tmp = join_cols(beta(set.A_), ones(1) * bound);
        iter_in = this->conjugate_gradient(tmp, XX, b, accuracy, max_iter) ;
      }
      beta(B) = theta * tmp.tail(1) ;
      beta(set.A_) = tmp.subvec(0,tmp.n_elem-2) ;
    }
    
    // Handling guys reaching the boundary (leaving the active set)
    double bound = penalty_.pen_norm(beta(B)) ; // current boundary
    uvec ind_toB = find(abs(beta(set.A_)) > bound) ;
    if (!ind_toB.is_empty()) {
      uvec toB = set.A_(ind_toB) ;
      set.del_vars(ind_toB) ;
      B = join_cols(B,toB);
      beta(B) = bound * sign(beta(B));
      theta = sign(beta(B));
    } else {
      success = true ;
    }
  }
  
  // Guys leaving the boundary after optimization (activation)
  grad = -data.XTy_ + XTX * beta ;
  uvec ind_toA  = find(theta == sign(grad(B)));
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

