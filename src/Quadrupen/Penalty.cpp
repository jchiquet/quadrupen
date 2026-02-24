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
vec Penalty<Norm::L1>::elt_norm(vec x, uvec pk) {
  return(arma::abs(x));
}

template<>
double Penalty<Norm::L1>::pen_norm(vec x, uvec pk) {
  return(sum(elt_norm(x)));
}

template<>
double Penalty<Norm::L1>::dual_norm(vec x, uvec pk) {
  return(max(elt_norm(x))) ;
}

template<>
vec Penalty<Norm::L1>::proximal(vec x, double lambda, uvec pk) {
  return(max(zeros(x.n_elem), 1-lambda/elt_norm(x)) % x);
}

// ______________________________________________________
// L2 NORM SQUARED A.K.A RIDGE
template<>
vec Penalty<Norm::RIDGE>::elt_norm(vec x, uvec pk) {
  return(arma::pow(x, 2));
}

template<>
double Penalty<Norm::RIDGE>::pen_norm(vec x, uvec pk) {
  return(sum(elt_norm(x)));
}

// Slight abuse: take the Fenchel conjuguate of L2^2
template<>
double Penalty<Norm::RIDGE>::dual_norm(vec x, uvec pk) {
  return(.25*sum(elt_norm(x))) ;
}

template<>
vec Penalty<Norm::RIDGE>::proximal(vec x, double lambda, uvec pk) {
  return(x / (1+2*lambda));
}

// ______________________________________________________
// LINF NORM A.K.A BOUNDED REGRESSION
template<>
vec Penalty<Norm::LINF>::elt_norm(vec x, uvec pk) {
  return(arma::abs(x));
}

template<>
double Penalty<Norm::LINF>::pen_norm(vec x, uvec pk) {
  return(max(elt_norm(x)));
}

template<>
double Penalty<Norm::LINF>::dual_norm(vec x, uvec pk) {
  return(sum(elt_norm(x))) ;
}

template<>
vec Penalty<Norm::LINF>::proximal(vec x, double lambda, uvec pk) {
  uword p = x.n_elem;
  vec u, proj;
  vec res = zeros<vec>(p);
  
  if ( as_scalar(sum(abs(x) / lambda)) >= 1) {
    
    // Reordering absolute values
    u = sort(abs(x), "descend");
    
    // values of the projected coordinate if non zero (dual problem)
    proj = (cumsum(u) - lambda)/linspace<vec>(1,p,p);
    
    // selecting nonnull entries (dual)
    uvec maxs = sort(find(u-proj>ZERO), "descend") ;
    double thresh = proj[maxs[0]];
    
    // solving primal problem
    // We keep the smallest values and threshold the common values to +- thresh
    for (uword k=0; k < p;k++) {
      if (fabs(x(k)) > ZERO) {
        if (x(k) > 0) {
          res(k) = fmin(fabs(x(k)),thresh);
        } else {
          res(k) = -fmin(fabs(x(k)),thresh);
        }
      }
    }
  }
  return(res);  
}

// ______________________________________________________
// L1/L2 NORM A.K.A GROUP-LASSO
template<>
vec Penalty<Norm::L1L2>::elt_norm(vec x, uvec pk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), 2);
    ind += pk(k);
  }
  
  return(res);
}

template<>
double Penalty<Norm::L1L2>::pen_norm(vec x, uvec pk) {
  return(sum(elt_norm(x, pk)));
}

template<>
double Penalty<Norm::L1L2>::dual_norm(vec x, uvec pk) {
  return(max(elt_norm(x, pk))) ;
}

template<>
vec Penalty<Norm::L1L2>::proximal(vec x, double lambda, uvec pk) {
  
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
vec Penalty<Norm::L1LINF>::elt_norm(vec x, uvec pk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), "inf");
    ind += pk(k);
  }
  
  return(res);
}

template<>
double Penalty<Norm::L1LINF>::pen_norm(vec x, uvec pk) {
  return(sum(elt_norm(x, pk)));
}

template<>
double Penalty<Norm::L1LINF>::dual_norm(vec x, uvec pk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), 1);
    ind += pk(k);
  }
  
  return(max(res)) ;
}

template<>
vec Penalty<Norm::L1LINF>::proximal(vec x, double lambda, uvec pk) {
  
  uword ind = 0 ;
  vec u, v, proj;
  vec res = zeros<vec>(sum(pk));
  
  for (uword k=0; k<pk.n_elem; k++) {
    v = x.subvec(ind,ind+pk(k)-1) ;
    uword p = v.n_elem ;
    
    // proximal l-inf
    if ( as_scalar(sum(abs(v) / lambda)) >= 1) {
      // Reordering absolute values
      u = sort(abs(v), "descend");
      
      // values of the projected coordinate if non zero (dual problem)
      proj = (cumsum(u) - lambda)/linspace<vec>(1,p,p);
      
      // selecting nonnull entries (dual)
      uvec maxs = sort(find(u-proj>ZERO), "descend") ;
      double thresh = proj[maxs[0]];
      
      // solving primal problem
      // We keep the smallest values and threshold the common values to +- thresh
      for (uword j=0; j<p ;j++) {
        if (fabs(v(j)) > ZERO) {
          v(j) = fmin(fabs(v(j)),thresh);
        } else {
          v(j) = -fmin(fabs(v(j)),thresh);
        }
      }
    }
    
    res.subvec(ind,ind+pk(k)-1) =  v ;
    ind += pk(k);
  }
  return(res);
}

// ______________________________________________________
// COOP(ERATIVE) NORM A.K.A COOPERATIVE-LASSO
template<>
vec Penalty<Norm::COOP>::elt_norm(vec x, uvec pk) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(max(zeros(pk(k)),x.subvec(ind, ind + pk(k) - 1)), 2) + norm(min(zeros(pk(k)),x.subvec(ind, ind + pk(k) - 1)),2);
    ind += pk(k);
  }
  
  return(res);
}

template<>
double Penalty<Norm::COOP>::pen_norm(vec x, uvec pk) {
  return(sum(elt_norm(x, pk)));
}

template<>
double Penalty<Norm::COOP>::dual_norm(vec x, uvec pk) {
  return(max(elt_norm(x, pk)));
}

template<>
vec Penalty<Norm::COOP>::proximal(vec x, double lambda, uvec pk) {
  
  vec res = zeros<vec>(x.n_elem);
  uword ind = 0 ;
  
  vec tmp = max(zeros(pk.n_elem), 1-lambda/elt_norm(x,pk)) ;
  
  for (uword k=0; k<pk.n_elem; k++) {
    res.subvec(ind, ind + pk(k) - 1) = tmp(k) * x.subvec(ind, ind + pk(k) - 1);
    ind += pk(k);
  }
  
  return(res);
}
