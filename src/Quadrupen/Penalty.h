#ifndef _quadrupen_PENALTY_H
#define _quadrupen_PENALTY_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>

#define ZERO 2e-16 // practical zero

enum class Norm {L1, LINF, RIDGE, L1L2, L1LINF, COOP, FUSED};

using namespace Rcpp;
using namespace arma;
using namespace std;

template <Norm norm>
class Penalty {
public: 
  
  Penalty() ;
  Penalty(SEXP PK) ;

  vec    elt_norm  (vec x) ;
  double pen_norm  (vec x) ;
  double dual_norm (vec x) ;
  vec proximal(vec x, double lambda) ;

  uvec pk ;
};

// Empty constructor
template <Norm norm>
Penalty<norm>::Penalty() {}

// Constructor requiring group information
template <Norm norm>
Penalty<norm>::Penalty(SEXP PK) : pk (as<uvec>(PK)){}

// ______________________________________________________
// L1 NORM A.K.A LASSO
template<>
inline vec Penalty<Norm::L1>::elt_norm(vec x) {
  return(arma::abs(x));
}

template<>
inline double Penalty<Norm::L1>::pen_norm(vec x) {
  return(sum(elt_norm(x)));
}

template<>
inline double Penalty<Norm::L1>::dual_norm(vec x) {
  return(max(elt_norm(x))) ;
}

template<>
inline vec Penalty<Norm::L1>::proximal(vec x, double lambda) {
  return(max(zeros(x.n_elem), 1-lambda/elt_norm(x)) % x);
}

// ______________________________________________________
// L2 NORM SQUARED A.K.A RIDGE
template<>
inline vec Penalty<Norm::RIDGE>::elt_norm(vec x) {
  return(arma::pow(x, 2));
}

template<>
inline double Penalty<Norm::RIDGE>::pen_norm(vec x) {
  return(sum(elt_norm(x)));
}

// Slight abuse: take the Fenchel conjuguate of L2^2
template<>
inline double Penalty<Norm::RIDGE>::dual_norm(vec x) {
  return(.25*sum(elt_norm(x))) ;
}

template<>
inline vec Penalty<Norm::RIDGE>::proximal(vec x, double lambda) {
  return(x / (1+2*lambda));
}

// ______________________________________________________
// LINF NORM A.K.A BOUNDED REGRESSION
template<>
inline vec Penalty<Norm::LINF>::elt_norm(vec x) {
  return(arma::abs(x));
}

template<>
inline double Penalty<Norm::LINF>::pen_norm(vec x) {
  return(max(elt_norm(x)));
}

template<>
inline double Penalty<Norm::LINF>::dual_norm(vec x) {
  return(sum(elt_norm(x))) ;
}

template<>
inline vec Penalty<Norm::LINF>::proximal(vec x, double lambda) {
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
inline vec Penalty<Norm::L1L2>::elt_norm(vec x) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), 2);
    ind += pk(k);
  }
  
  return(res);
}

template<>
inline double Penalty<Norm::L1L2>::pen_norm(vec x) {
  return(sum(elt_norm(x)));
}

template<>
inline double Penalty<Norm::L1L2>::dual_norm(vec x) {
  return(max(elt_norm(x))) ;
}

template<>
inline vec Penalty<Norm::L1L2>::proximal(vec x, double lambda) {
  
  vec res = zeros<vec>(x.n_elem);
  uword ind = 0 ;
  
  vec tmp = max(zeros(pk.n_elem), 1-lambda/elt_norm(x)) ;
  
  for (uword k=0; k<pk.n_elem; k++) {
    res.subvec(ind, ind + pk(k) - 1) = tmp(k) * x.subvec(ind, ind + pk(k) - 1);
    ind += pk(k);
  }
  
  return(res);
}

// ______________________________________________________
// L1/LINF NORM A.K.A GROUP-LASSO type 2
template<>
inline vec Penalty<Norm::L1LINF>::elt_norm(vec x) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), "inf");
    ind += pk(k);
  }
  
  return(res);
}

template<>
inline double Penalty<Norm::L1LINF>::pen_norm(vec x) {
  return(sum(elt_norm(x)));
}

template<>
inline double Penalty<Norm::L1LINF>::dual_norm(vec x) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), 1);
    ind += pk(k);
  }
  
  return(max(res)) ;
}

template<>
inline vec Penalty<Norm::L1LINF>::proximal(vec x, double lambda) {
  
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
inline vec Penalty<Norm::COOP>::elt_norm(vec x) {
  
  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  uword ind = 0 ; // index to go through the groups
  
  for (uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(max(zeros(pk(k)),x.subvec(ind, ind + pk(k) - 1)), 2) + norm(min(zeros(pk(k)),x.subvec(ind, ind + pk(k) - 1)),2);
    ind += pk(k);
  }
  
  return(res);
}

template<>
inline double Penalty<Norm::COOP>::pen_norm(vec x) {
  return(sum(elt_norm(x)));
}

template<>
inline double Penalty<Norm::COOP>::dual_norm(vec x) {
  return(max(elt_norm(x)));
}

template<>
inline vec Penalty<Norm::COOP>::proximal(vec x, double lambda) {
  
  vec res = zeros<vec>(x.n_elem);
  uword ind = 0 ;
  
  vec tmp = max(zeros(pk.n_elem), 1-lambda/elt_norm(x)) ;
  
  for (uword k=0; k<pk.n_elem; k++) {
    res.subvec(ind, ind + pk(k) - 1) = tmp(k) * x.subvec(ind, ind + pk(k) - 1);
    ind += pk(k);
  }
  
  return(res);
}

#endif
