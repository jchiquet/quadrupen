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

#endif
