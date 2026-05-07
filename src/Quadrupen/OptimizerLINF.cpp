/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "OptimizerLINF.h"

using namespace Rcpp;
using namespace arma;

OptimizerLINF::OptimizerLINF(
  SimplePenalty<SimpleNorm::LINF>& penalty, const List& control) : 
  Optimizer(control) 
{penalty_ = penalty ;}

mat OptimizerLINF::updateCholeskyFromExisting(const mat& R, const vec& b) {
  uword p = R.n_cols + 1;
  mat R_ ;
  colvec rp  = zeros<colvec>(p,1);
  rp.subvec(0,p-2) = arma::solve (trimatl(strans(R)), b.subvec(0,p-2));
  rp(p-1) = sqrt(b(p-1) - dot(rp,rp));
  R_ = join_rows( join_cols(R_, zeros<mat>(1,p-1)) , rp);
  return(R_) ;
}

uword OptimizerLINF::quadratic_breg(
  vec& beta,
  vec &grad,
  const double& lambda,
  const vec& weights,
  RegressionData<mat> &data,
  uvec& active,
  const double& accuracy,
  const uword& max_iter) {
  
  uvec B = regspace<uvec>(0,beta.n_elem-1) ; B.shed_rows(active) ;
  grad = -data.XTy_ + data.XTX_ * beta ;
  vec theta = -sign(grad(B)) ; // sign of the guys on the boundary
  
  uword iter = 0, iter_in= 0 ; // count the number of systems solved
  bool success = false ;
  while ((!success) && (iter < 3)) {
    
    iter++;
    
    // SOLVE THE QUADRATIC PROBLEM
    vec XX_B = data.XTX_.cols(B) * theta;
    if (active.is_empty()) {
      double b  = dot(theta, data.XTy_(B)) - lambda * mean(weights(B));
      beta(B) = theta * (b/sum(theta % XX_B(B),0)) ;
    } else {
      // Constructing the system (KKT)
      vec tmp = join_cols(XX_B(active),sum(theta % XX_B(B),0)) ;
      mat XX = join_rows(join_cols(data.XTX_(active, active), trans(XX_B(active))), tmp);
      vec b = data.XTy_(active) ; b.resize(b.n_elem + 1) ;
      b.tail(1) = dot(theta, data.XTy_(B)) - lambda * mean(weights(B));
      vec pen_arg = beta(B);
      double bound = max(abs(pen_arg)) ;
      tmp = join_cols(beta(active), ones(1) * bound);
      iter_in = this->conjugate_gradient(tmp, XX, b, accuracy, max_iter) ;
      beta(B) = theta * tmp.tail(1) ;
      beta(active) = tmp.subvec(0,tmp.n_elem-2) ;
    }
    
    // Handling guys reaching the boundary (leaving the active set)
    vec pen_arg = beta(B);
    double bound = max(abs(pen_arg)) ; // current boundary
    uvec ind_toB = find(abs(beta(active)) > bound) ;
    if (!ind_toB.is_empty()) {
      uvec toB = active(ind_toB) ;
      active.shed_rows(ind_toB) ;
      B = join_cols(B, toB);
      beta(B) = bound * sign(beta(B));
      theta = sign(beta(B));
    } else {
      success = true ;
    }
  }
  
  // Guys leaving the boundary after optimization (activation)
  grad = -data.XTy_ + data.XTX_ * beta ;
  uvec ind_toA  = find(theta == sign(grad(B)));
  if (!ind_toA.is_empty()) {
    uvec toA = B(ind_toA) ;
    active.resize(active.n_elem + toA.n_elem);
    active.tail(toA.n_elem) = toA;
    B.shed_rows(ind_toA) ;
    if (B.is_empty()) {
      throw std::runtime_error("Every variable left the boudary. Try more regularization.");
    }
  }
  return(iter_in) ;
  
}
