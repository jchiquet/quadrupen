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

  vec x_st = x;
  if (alpha_ > 0) 
    x_st = sign(x_st) % max(abs(x_st) - alpha_, zeros<vec>(x_st.n_elem)) ;

  vec res = zeros<vec>(pk.n_elem);
  uword ind = 0;
  for (uword k = 0; k < pk.n_elem; k++) {
    uword size = pk(k);
    res(k) = norm(x_st.rows(ind, ind + size - 1), 2);
    ind += size;
  }

  return(res / (wk  * (1 - alpha_))) ;
}

template<>
vec GroupPenalty<GroupNorm::L1L2>::proximal(const vec& x, double lambda, const uvec& pk, const vec& wk) {

  vec res = x;
  
  double l1_threshold = lambda * alpha_;
  double l2_factor = lambda * (1.0 - alpha_);
  
  if (l1_threshold > 0) {
    res = sign(x) % max(abs(x) - l1_threshold, zeros<vec>(x.n_elem));
  }
  
  uword ind = 0;
  for (uword k = 0; k < pk.n_elem; k++) {
    uword group_size = pk(k);
    auto group_slice = res.subvec(ind, ind + group_size - 1);
    
    double group_norm = norm(group_slice, 2);
    double threshold_g = l2_factor * wk(k);
    if (group_norm > threshold_g && group_norm > 0) {
      double shrink = 1.0 - (threshold_g / group_norm);
      group_slice *= shrink;
    } else {
      group_slice.zeros();
    }
    ind += group_size;
  }
  
  return res;
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
  
  vec x_st = sign(x) % max(abs(x) - alpha_, zeros<vec>(x.n_elem)) ;

  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = sum(abs(x_st.subvec(ind, ind + pk(k) - 1))) ;
    ind += pk(k);
  }
  
  return(res / (wk * (1-alpha_))) ;
  
}

template<>
vec GroupPenalty<GroupNorm::L1LINF>::proximal(const vec& x, double lambda, const uvec& pk, const vec& wk) {
  
  vec res ;
  if (alpha_ > 0) {
    res = sign(x) % max(abs(x) - alpha_, zeros<vec>(x.n_elem)) ;
  } else {
    res = x ;
  }
  
  uword ind = 0 ;
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
  
  vec x_st = sign(x) % max(abs(x) - alpha_, zeros<vec>(x.n_elem)) ;
  
  vec  res = zeros<vec> (pk.n_elem) ;
  uword ind = 0 ; // index to go through the groups
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
    res = sign(x) % max(abs(x) - alpha_, zeros<vec>(x.n_elem)) ;
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

