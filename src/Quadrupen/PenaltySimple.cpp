/*
  * Author: Julien CHIQUET
*         MIA Paris-Saclay
*/
  
#include "PenaltySimple.h"

using namespace Rcpp;
using namespace arma;

// ______________________________________________________
// L1 NORM A.K.A LASSO
template<>
vec SimplePenalty<SimpleNorm::L1>::elt_norm(const vec& x) {
  return(arma::abs(x));
}

template<>
vec SimplePenalty<SimpleNorm::L1>::elt_dual_norm(const vec& x) {
  return(arma::abs(x));
}

template<>
double SimplePenalty<SimpleNorm::L1>::pen_norm(const vec& x, const vec& w) {
  return(accu(w % elt_norm(x)));
}

template<>
double SimplePenalty<SimpleNorm::L1>::dual_norm(const vec& x, const vec& w) {
  return(max(elt_dual_norm(x) / w)) ;
}

template<>
vec SimplePenalty<SimpleNorm::L1>::proximal(const vec& x, double lambda, const vec& w) {
  return(x % clamp(1 - lambda * w/elt_norm(x), 0, datum::inf));
}

// ______________________________________________________
// L2 NORM SQUARED A.K.A RIDGE
template<>
vec SimplePenalty<SimpleNorm::RIDGE>::elt_norm(const vec& x) {
  return(arma::pow(x, 2));
}

template<>
vec SimplePenalty<SimpleNorm::RIDGE>::elt_dual_norm(const vec& x) {
  return(arma::pow(x, 2));
}

template<>
double SimplePenalty<SimpleNorm::RIDGE>::pen_norm(const vec& x, const vec& w) {
  return(sum(w % elt_norm(x)));
}

// Slight abuse: take the Fenchel conjugate of L2^2
template<>
double SimplePenalty<SimpleNorm::RIDGE>::dual_norm(const vec& x, const vec& w) {
  return(.25*sum(elt_norm(x) / w)) ;
}

template<>
vec SimplePenalty<SimpleNorm::RIDGE>::proximal(const vec& x, double lambda, const vec& w) {
  return(x / (1+2*lambda*w));
}

// ______________________________________________________
// LINF NORM A.K.A BOUNDED REGRESSION
template<>
vec SimplePenalty<SimpleNorm::LINF>::elt_norm(const vec& x) {
  return(arma::abs(x));
}

template<>
double SimplePenalty<SimpleNorm::LINF>::pen_norm(const vec& x, const vec& w) {
  return(max(elt_norm(x) % w));
}

template<>
double SimplePenalty<SimpleNorm::LINF>::dual_norm(const vec& x, const vec& w) {
  return(sum(elt_norm(x) / w)) ;
}

template<>
vec SimplePenalty<SimpleNorm::LINF>::proximal(const vec& x, double lambda, const vec& w) {
  uword p = x.n_elem;
  vec res = zeros<vec>(p);

  // Project onto the l1 ball
  if (accu(w % abs(x)) > lambda) {
      
    // Reordering absolute values
    vec u = sort(abs(x) / w, "descend");

    // values of the projected coordinate if non zero (dual problem)
    vec proj = (cumsum(u) - lambda)/linspace<vec>(1,p,p);    
    
    // find critical index
    uword i = max(find(u - proj >= 0)) ;
    
    res = sign(x) % min(abs(x) % w, proj(i) * ones(x.n_elem) );
  }
  return(res);  
}

