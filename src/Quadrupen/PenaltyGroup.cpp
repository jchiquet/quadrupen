/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "PenaltyGroup.h"

using namespace Rcpp;
using namespace arma;

// ______________________________________________________
// L1/L2 NORM A.K.A GROUP-LASSO
template<>
vec GroupPenalty<GroupNorm::L1L2>::elt_norm(const vec& x, const uvec& pk, const vec& wk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    if (alpha_ > 0) res(k) += alpha_ * sum(abs(x.subvec(ind, ind + pk(k) - 1))) ;
    res(k) += (1 - alpha_) * wk(k) * norm(x.subvec(ind, ind + pk(k) - 1), 2) ;
    ind += pk(k);
  }
  
  return(res);
}

template<>
vec GroupPenalty<GroupNorm::L1L2>::elt_dual_norm  (const vec& x, const uvec& pk, const vec& wk) {
  
  vec x_st = x % clamp(1 - alpha_ / abs(x), 0, datum::inf) ;
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) += norm(x_st.subvec(ind, ind + pk(k) - 1), 2);
    ind += pk(k);
  }
  
  return(res / (wk  * (1 - alpha_))) ;
}

template<>
vec GroupPenalty<GroupNorm::L1L2>::proximal(const vec& x, double lambda, const uvec& pk, const vec& wk) {
  
  vec res ;
  
  if (alpha_ > 0) {
    res = x % clamp(1 - lambda * alpha_ / abs(x), 0, datum::inf);
  } else {
    res = x ;
  }
  
  uword ind = 0 ;
  for (uword k=0; k<pk.n_elem; k++) {
    double shrink = fmax(0, 1 - lambda * (1-alpha_) * wk(k)/norm(res.subvec(ind, ind + pk(k) - 1), 2)) ;
    res.subvec(ind, ind + pk(k) - 1) *= shrink ;
    ind += pk(k);
  }
  
  return(res);
}

// ______________________________________________________
// L1/LINF NORM A.K.A GROUP-LASSO type 2
template<>
vec GroupPenalty<GroupNorm::L1LINF>::elt_norm(const vec& x, const uvec& pk, const vec& wk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    if (alpha_ > 0) res(k) += alpha_ * norm(x.subvec(ind, ind + pk(k) - 1), 1);
    res(k) += (1 - alpha_) * wk(k) * max(abs(x.subvec(ind, ind + pk(k) - 1)));
    ind += pk(k);
  }
  
  return(res);
}

template<>
vec GroupPenalty<GroupNorm::L1LINF>::elt_dual_norm  (const vec& x, const uvec& pk, const vec& wk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  vec x_st = x % clamp(1 - alpha_ / abs(x), 0, datum::inf) ;
  
  uword ind = 0 ; // index to go through the groups
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = sum(abs(x_st.subvec(ind, ind + pk(k) - 1))) ;
    ind += pk(k);
  }
  
  return(res / (wk * (1-alpha_))) ;
  
}

template<>
vec GroupPenalty<GroupNorm::L1LINF>::proximal(const vec& x, double lambda, const uvec& pk, const vec& wk) {
  
  uword ind = 0 ;
  vec res = x ;
  
  if (alpha_ > 0) {
    res = res % clamp(1 - alpha_ / abs(res), 0, datum::inf);
  }
  
  for (uword k=0; k<pk.n_elem; k++) {
    uword p = pk(k);
    vec res_g = res.subvec(ind,ind+p-1) ;
    
    if ( accu(abs(res_g)) > lambda * wk(k)) {
      // Project onto the l1 ball
      
      // Reordering absolute values
      vec u = sort(abs(res_g), "descend");
      
      // values of the projected coordinate if non zero (dual problem)
      vec proj = (cumsum(u) - lambda * wk(k))/linspace<vec>(1,p,p);
      
      // find critical index
      uword i = max(find(u - proj >= 0)) ;
      
      // res = max(zeros(x.n_elem), abs(x) - proj[i]) ;
      res_g = sign(res_g) % min(abs(res_g), proj(i) * ones(p) );
    }
    
    res.subvec(ind,ind+p-1) = (1 - alpha_) * res_g ;
    ind += p;
  }
  return(res);
}

// ______________________________________________________
// COOP(ERATIVE) NORM A.K.A COOPERATIVE-LASSO
template<>
vec GroupPenalty<GroupNorm::COOP>::elt_norm(const vec& x, const uvec& pk, const vec& wk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = wk(k) * (
      norm(max(zeros(pk(k)),   x.subvec(ind, ind + pk(k) - 1)), 2) + 
        norm(max(zeros(pk(k)), - x.subvec(ind, ind + pk(k) - 1)), 2)
    );
    ind += pk(k);
  }
  
  return(res);
}

template<>
vec GroupPenalty<GroupNorm::COOP>::elt_dual_norm(const vec& x, const uvec& pk, const vec& wk) {
  
  vec  res = zeros<vec> (pk.n_elem) ;
  uword ind = 0 ; // index to go through the groups
  
  vec x_st = x % clamp(1 - alpha_ / abs(x), 0, datum::inf) ;
  
  for (uword k=0; k<pk.n_elem; k++) {
    double norm_pos = norm(max(zeros(pk(k)),   x_st.subvec(ind, ind + pk(k) - 1)), 2) ;
    double norm_neg = norm(max(zeros(pk(k)), - x_st.subvec(ind, ind + pk(k) - 1)), 2) ;
    res(k) = (norm_pos > norm_neg) ? norm_pos : norm_neg ;
    res(k) /= wk(k) * (1-alpha_) ;
    ind += pk(k);
  }
  
  return(res);
}

template<>
vec GroupPenalty<GroupNorm::COOP>::proximal(const vec& x, double lambda, const uvec& pk, const vec& wk) {
  
  vec res ;
  
  if (alpha_ > 0) {
    res = x % clamp(1 - lambda * alpha_ / abs(x), 0, datum::inf);
  } else {
    res = x ;
  }
  
  uword ind = 0;
  
  for (uword k=0; k<pk.n_elem; k++) {
    double shrink_pos = fmax(0, 1 - lambda * (1-alpha_) * wk(k)/norm(max(zeros(pk(k)),   res.subvec(ind, ind + pk(k) - 1)), 2)) ;
    double shrink_neg = fmax(0, 1 - lambda * (1-alpha_) * wk(k)/norm(max(zeros(pk(k)), - res.subvec(ind, ind + pk(k) - 1)), 2));
    
    for (uword j=ind; j<(ind+pk(k)); j++) {
      if(res[j] > 0) {
        res[j] *= shrink_pos ;
      } else {
        res[j] *= shrink_neg ;	
      }
    }
    ind += pk(k);
  }
  
  return(res);
}

