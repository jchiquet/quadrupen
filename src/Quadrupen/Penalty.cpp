/*
  * Author: Julien CHIQUET
*         MIA Paris-Saclay
*/
  
#include "Penalty.h"

using namespace Rcpp;
using namespace arma;

// ______________________________________________________
// L1 NORM A.K.A LASSO
template<>
vec SimplePenalty<SimpleNorm::L1>::elt_norm(vec x) {
  return(arma::abs(x));
}

template<>
double SimplePenalty<SimpleNorm::L1>::pen_norm(vec x) {
  return(accu(elt_norm(x)));
}

template<>
double SimplePenalty<SimpleNorm::L1>::dual_norm(vec x) {
  return(max(elt_norm(x))) ;
}

template<>
vec SimplePenalty<SimpleNorm::L1>::proximal(vec x, double lambda) {
  return(max(zeros(x.n_elem), 1-lambda/elt_norm(x)) % x);
}

// ______________________________________________________
// L2 NORM SQUARED A.K.A RIDGE
template<>
vec SimplePenalty<SimpleNorm::RIDGE>::elt_norm(vec x) {
  return(arma::pow(x, 2));
}

template<>
double SimplePenalty<SimpleNorm::RIDGE>::pen_norm(vec x) {
  return(sum(elt_norm(x)));
}

// Slight abuse: take the Fenchel conjuguate of L2^2
template<>
double SimplePenalty<SimpleNorm::RIDGE>::dual_norm(vec x) {
  return(.25*sum(elt_norm(x))) ;
}

template<>
vec SimplePenalty<SimpleNorm::RIDGE>::proximal(vec x, double lambda) {
  return(x / (1+2*lambda));
}

// ______________________________________________________
// LINF NORM A.K.A BOUNDED REGRESSION
template<>
vec SimplePenalty<SimpleNorm::LINF>::elt_norm(vec x) {
  return(arma::abs(x));
}

template<>
double SimplePenalty<SimpleNorm::LINF>::pen_norm(vec x) {
  return(max(elt_norm(x)));
}

template<>
double SimplePenalty<SimpleNorm::LINF>::dual_norm(vec x) {
  return(sum(elt_norm(x))) ;
}

template<>
vec SimplePenalty<SimpleNorm::LINF>::proximal(vec x, double lambda) {
  uword p = x.n_elem;
  vec res = zeros<vec>(p);
  
  // Project onto the l1 ball
  if (accu(abs(x)) > lambda) {
    // Project onto the l1 ball

    // Reordering absolute values
    vec u = sort(abs(x), "descend");
    
    // values of the projected coordinate if non zero (dual problem)
    vec proj = (cumsum(u) - lambda)/linspace<vec>(1,p,p);
    
    // find critical index
    uword i = max(find(u - proj >= 0)) ;
    
    // res = max(zeros(x.n_elem), abs(x) - proj[i]) ;
    res = sign(x) % min(abs(x), proj(i) * ones(x.n_elem) );
  }
  return(res);  
}

// ______________________________________________________
// L1/L2 NORM A.K.A GROUP-LASSO
template<>
vec MixedPenalty<MixedNorm::L1L2>::elt_norm(vec x, uvec pk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), 2);
    ind += pk(k);
  }
  
  return(res);
}

template<>
double MixedPenalty<MixedNorm::L1L2>::pen_norm(vec x, uvec pk, vec wk) {
  return(sum(wk % elt_norm(x, pk)));
}

template<>
vec MixedPenalty<MixedNorm::L1L2>::elt_dual_norm  (vec x, uvec pk) {
  return(elt_norm(x, pk)) ;
}

template<>
double MixedPenalty<MixedNorm::L1L2>::dual_norm(vec x, uvec pk, vec wk) {
  return(max(elt_dual_norm(x, pk) / wk)) ;
}

template<>
vec MixedPenalty<MixedNorm::L1L2>::proximal(vec x, vec lambda, uvec pk) {
  
  vec res = zeros<vec>(x.n_elem);
  uword ind = 0 ;
  
  vec tmp = max(zeros(pk.n_elem), 1-lambda/elt_norm(x,pk)) ;
  
  for (uword k=0; k<pk.n_elem; k++) {
    res.subvec(ind, ind + pk(k) - 1) = tmp(k) * x.subvec(ind, ind + pk(k) - 1);
    ind += pk(k);
  }
  
  return(res);
}

// ______________________________________________________
// L1/LINF NORM A.K.A GROUP-LASSO type 2
template<>
vec MixedPenalty<MixedNorm::L1LINF>::elt_norm(vec x, uvec pk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = max(abs(x.subvec(ind, ind + pk(k) - 1)));
    ind += pk(k);
  }
  
  return(res);
}

template<>
double MixedPenalty<MixedNorm::L1LINF>::pen_norm(vec x, uvec pk, vec wk) {
  return(sum(wk % elt_norm(x, pk)));
}

template<>
vec MixedPenalty<MixedNorm::L1LINF>::elt_dual_norm  (vec x, uvec pk) {

  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = sum(abs(x.subvec(ind, ind + pk(k) - 1)));
    ind += pk(k);
  }
  
  return(res) ;
  
}

template<>
double MixedPenalty<MixedNorm::L1LINF>::dual_norm(vec x, uvec pk, vec wk) {
  return(max(elt_dual_norm(x, pk) / wk)) ;
}

template<>
vec MixedPenalty<MixedNorm::L1LINF>::proximal(vec x, vec lambda, uvec pk) {

  uword ind = 0 ;
  vec res = zeros<vec>(x.n_elem);
  
  for (uword k=0; k<pk.n_elem; k++) {
    uword p = pk(k);
    vec x_g = x.subvec(ind,ind+p-1) ;

    if ( accu(abs(x_g)) > lambda(k)) {
      // Project onto the l1 ball
      
      // Reordering absolute values
      vec u = sort(abs(x_g), "descend");
      
      // values of the projected coordinate if non zero (dual problem)
      vec proj = (cumsum(u) - lambda(k))/linspace<vec>(1,p,p);
      
      // find critical index
      uword i = max(find(u - proj >= 0)) ;
      
      // res = max(zeros(x.n_elem), abs(x) - proj[i]) ;
      x_g = sign(x_g) % min(abs(x_g), proj(i) * ones(p) );
    }

    res.subvec(ind,ind+p-1) = x_g ;
    ind += p;
  }
  return(res);
}

// ______________________________________________________
// COOP(ERATIVE) NORM A.K.A COOPERATIVE-LASSO
template<>
vec MixedPenalty<MixedNorm::COOP>::elt_norm(vec x, uvec pk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = 
      norm(max(zeros(pk(k)),   x.subvec(ind, ind + pk(k) - 1)), 2) + 
      norm(max(zeros(pk(k)), - x.subvec(ind, ind + pk(k) - 1)), 2);
    ind += pk(k);
  }
  
  return(res);
}

template<>
double MixedPenalty<MixedNorm::COOP>::pen_norm(vec x, uvec pk, vec wk) {
  return(sum(wk % elt_norm(x, pk)));
}

template<>
vec MixedPenalty<MixedNorm::COOP>::elt_dual_norm(vec x, uvec pk) {

  vec  res = zeros<vec> (pk.n_elem) ;
  uword ind = 0 ; // index to go through the groups

  for (uword k=0; k<pk.n_elem; k++) {
    double norm_pos = norm(max(zeros(pk(k)),   x.subvec(ind, ind + pk(k) - 1)), 2) ;
    double norm_neg = norm(max(zeros(pk(k)), - x.subvec(ind, ind + pk(k) - 1)), 2) ;
    res(k) = (norm_pos > norm_neg) ? norm_pos : norm_neg ;
    ind += pk(k);
  }

  return(res);
}
  
template<>
double MixedPenalty<MixedNorm::COOP>::dual_norm(vec x, uvec pk, vec wk) {
  return(max(elt_dual_norm(x, pk) / wk));
}

template<>
vec MixedPenalty<MixedNorm::COOP>::proximal(vec x, vec lambda, uvec pk) {
  
  vec res = x ;
  uword ind = 0;
  
  for (uword k=0; k<pk.n_elem; k++) {
    double shrink_pos = fmax(0, 1 - lambda(k)/norm(max(zeros(pk(k)),   res.subvec(ind, ind + pk(k) - 1)), 2)) ;
    double shrink_neg = fmax(0, 1 - lambda(k)/norm(max(zeros(pk(k)), - res.subvec(ind, ind + pk(k) - 1)), 2));
    
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
