/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "PenaltyDense.h"

using namespace Rcpp;
using namespace arma;

// ______________________________________________________
// LINF NORM A.K.A BOUNDED REGRESSION
template<>
vec DensePenalty<DenseNorm::LINF>::elt_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) % w);
}

template<>
vec DensePenalty<DenseNorm::LINF>::elt_dual_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) / w);
}

template<>
double DensePenalty<DenseNorm::LINF>::pen_norm(const vec& x, const vec& w, double lambda) {
  return(max(elt_norm(x, w)));
}

template<>
double DensePenalty<DenseNorm::LINF>::dual_norm(const vec& x, const vec& w, double lambda) {
  return(sum(elt_dual_norm(x,w))) ;
}

template<>
vec DensePenalty<DenseNorm::LINF>::proximal(const vec& x, double lambda, const vec& w) {
  uword p = x.n_elem;
  vec abs_x_w = arma::abs(x) / w;
  
  if (accu(abs_x_w) <= lambda) {
    return zeros<vec>(p);
  }
  
  // Project onto the l1 ball
  
  // Reordering absolute values
  vec u = sort(abs_x_w, "descend");
  
  // values of the projected coordinate if non zero (dual problem)
  // vec proj = (cumsum(u) - lambda)/linspace<vec>(1,p,p);
  // 
  // find critical index
  // uword i = max(find(u - proj >= 0)) ;
  
  double cumsum_u_w = 0.0;
  double cumsum_w_sq = 0.0;
  double rho = 0.0;
  
  for (uword j = 0; j < p; ++j) {
    
    cumsum_u_w += u(j) * (w(j) * w(j)); //
    cumsum_w_sq += w(j) * w(j);
    
    double t = (cumsum_u_w - lambda) / cumsum_w_sq;
    
    if (j < p - 1 && u(j + 1) <= t) {
      rho = t;
      break;
    }
    if (j == p - 1) rho = t;
  }
  
  vec res = sign(x) % min(abs(x), rho * w );
  
  return(res);
}

// ______________________________________________________
// L2 NORM SQUARED A.K.A RIDGE
template<>
vec DensePenalty<DenseNorm::L2>::elt_norm(const vec& x, const vec& w, double lambda) {
  return(arma::pow(x, 2) % w);
}

template<>
vec DensePenalty<DenseNorm::L2>::elt_dual_norm(const vec& x, const vec& w, double lambda) {
  return(arma::pow(x, 2) / w);
}

template<>
double DensePenalty<DenseNorm::L2>::pen_norm(const vec& x, const vec& w, double lambda) {
  return(sum(elt_norm(x, w)));
}

// Slight abuse: take the Fenchel conjugate of L2^2
template<>
double DensePenalty<DenseNorm::L2>::dual_norm(const vec& x, const vec& w, double lambda) {
  return(.25*sum(elt_norm(x, w))) ;
}

template<>
vec DensePenalty<DenseNorm::L2>::proximal(const vec& x, double lambda, const vec& w) {
  return(x / (1+2*lambda*w));
}
