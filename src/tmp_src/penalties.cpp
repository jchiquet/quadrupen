/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#include "penalties.hpp"

using namespace Rcpp;
using namespace arma;

// Empty constructor for basic penalties
PENALTY::PENALTY() {}

// Constructor requiring group information
PENALTY::PENALTY(SEXP PK) : pk (as<uvec>(PK)){}

void PENALTY::setPenalty(std::string penName){

  // link public member functions to l1 norm functions
  if (penName ==  "l1") {
    _current_elt_norm_ptr  = &PENALTY::elt_norm_L1;
    _current_pen_norm_ptr  = &PENALTY::pen_norm_L1;
    _current_dual_norm_ptr = &PENALTY::dual_norm_L1;
    _current_proximal_ptr  = &PENALTY::proximal_L1;

  // link public member functions to linf norm functions
  } else if (penName == "linf") {
    _current_elt_norm_ptr  = &PENALTY::elt_norm_LINF;
    _current_pen_norm_ptr  = &PENALTY::pen_norm_LINF;
    _current_dual_norm_ptr = &PENALTY::dual_norm_LINF;
    _current_proximal_ptr  = &PENALTY::proximal_LINF;

  // link public member functions to l1/l2 norm functions
  } else if (penName ==  "l1l2") {
    _current_elt_norm_ptr  = &PENALTY::elt_norm_L1L2;
    _current_pen_norm_ptr  = &PENALTY::pen_norm_L1L2;
    _current_dual_norm_ptr = &PENALTY::dual_norm_L1L2;
    _current_proximal_ptr  = &PENALTY::proximal_L1L2;

  // link public member functions to l1/linf norm functions
  } else if (penName == "l1linf") {
    _current_elt_norm_ptr  = &PENALTY::elt_norm_L1LINF;
    _current_pen_norm_ptr  = &PENALTY::pen_norm_L1LINF;
    _current_dual_norm_ptr = &PENALTY::dual_norm_L1LINF;
    _current_proximal_ptr  = &PENALTY::proximal_L1LINF;

  // link public member functions to coop norm functions
  } else if (penName == "coop") {
    _current_elt_norm_ptr  = &PENALTY::elt_norm_COOP;
    _current_pen_norm_ptr  = &PENALTY::pen_norm_COOP;
    _current_dual_norm_ptr = &PENALTY::dual_norm_COOP;
    _current_proximal_ptr  = &PENALTY::proximal_COOP;
  }
}

vec    PENALTY::elt_norm  (vec x) {return((this->*_current_elt_norm_ptr)(x));}
double PENALTY::pen_norm  (vec x) {return((this->*_current_pen_norm_ptr)(x));}
double PENALTY::dual_norm (vec x) {return((this->*_current_dual_norm_ptr)(x));}
vec    PENALTY::proximal  (vec x, double lambda) {return((this->*_current_proximal_ptr)(x,lambda));}

// ______________________________________________________
// L1 NORM A.K.A LASSO
vec PENALTY::elt_norm_L1(vec x) {
  return(arma::abs(x));
}

double PENALTY::pen_norm_L1(vec x) {
  return(sum(elt_norm(x)));
}

double PENALTY::dual_norm_L1(vec x) {
  return(max(elt_norm(x))) ;
}

vec PENALTY::proximal_L1(vec x, double lambda) {
  return(max(zeros(x.n_elem), 1-lambda/elt_norm(x)) % x);
}

// ______________________________________________________
// LINF NORM A.K.A BOUNDED REGRESSION
vec PENALTY::elt_norm_LINF(vec x) {
  return(arma::abs(x));
}

double PENALTY::pen_norm_LINF(vec x) {
  return(max(elt_norm(x)));
}

double PENALTY::dual_norm_LINF(vec x) {
  return(sum(elt_norm(x))) ;
}

vec PENALTY::proximal_LINF(vec x, double lambda) {
  int p = x.n_elem;
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
    for (int k=0; k < p;k++) {
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
vec PENALTY::elt_norm_L1L2(vec x) {

  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  int ind = 0 ; // index to go through the groups

  for (int k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), 2);
    ind += pk(k);
  }

  return(res);
}

double PENALTY::pen_norm_L1L2(vec x) {
  return(sum(elt_norm(x)));
}

double PENALTY::dual_norm_L1L2(vec x) {
  return(max(elt_norm(x))) ;
}

vec PENALTY::proximal_L1L2(vec x, double lambda) {

  vec res = zeros<vec>(x.n_elem);
  int k,ind = 0 ;
  
  vec tmp = max(zeros(pk.n_elem), 1-lambda/elt_norm(x)) ;

  for (k=0; k<pk.n_elem; k++) {
    res.subvec(ind, ind + pk(k) - 1) = tmp(k) * x.subvec(ind, ind + pk(k) - 1);
    ind += pk(k);
  }

  return(res);
}

// ______________________________________________________
// L1/LINF NORM A.K.A GROUP-LASSO type 2
vec PENALTY::elt_norm_L1LINF(vec x) {

  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  int ind = 0 ; // index to go through the groups

  for (int k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), "inf");
    ind += pk(k);
  }
  
  return(res);
}

double PENALTY::pen_norm_L1LINF(vec x) {
  return(sum(elt_norm(x)));
}

double PENALTY::dual_norm_L1LINF(vec x) {

  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  int ind = 0 ; // index to go through the groups

  for (int k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), 1);
    ind += pk(k);
  }
  
  return(max(res)) ;
}

vec PENALTY::proximal_L1LINF(vec x, double lambda) {

  uword ind = 0, p ;
  vec u, v, proj;
  vec res = zeros<vec>(sum(pk));

  for (int k=0; k<pk.n_elem; k++) {
    v = x.subvec(ind,ind+pk(k)-1) ;
    p = v.n_elem ;
    
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
      for (int j=0; j<p ;j++) {
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
vec PENALTY::elt_norm_COOP(vec x) {

  vec  res = zeros<vec> (pk.n_elem) ; // output with group norms
  int ind = 0 ; // index to go through the groups

  for (int k=0; k<pk.n_elem; k++) {
    res(k) = norm(max(zeros(pk(k)),x.subvec(ind, ind + pk(k) - 1)), 2) + norm(min(zeros(pk(k)),x.subvec(ind, ind + pk(k) - 1)),2);
    ind += pk(k);
  }

  return(res);
}

double PENALTY::pen_norm_COOP(vec x) {
  return(sum(elt_norm(x)));
}

double PENALTY::dual_norm_COOP(vec x) {
  return(max(elt_norm(x)));
}

vec PENALTY::proximal_COOP(vec x, double lambda) {

  vec res = zeros<vec>(x.n_elem);
  int k,ind = 0 ;
  
  vec tmp = max(zeros(pk.n_elem), 1-lambda/elt_norm(x)) ;

  for (k=0; k<pk.n_elem; k++) {
    res.subvec(ind, ind + pk(k) - 1) = tmp(k) * x.subvec(ind, ind + pk(k) - 1);
    ind += pk(k);
  }

  return(res);

}

