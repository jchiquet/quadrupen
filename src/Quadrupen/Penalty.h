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

enum class Norm {L1, LINF, RIDGE, L1L2, L1LINF, COOP};

using namespace Rcpp;
using namespace arma;
using namespace std;

template <Norm norm>
class Penalty {
public: 
  
  Penalty() {} ;

  double pen_norm  (vec x, uvec pk = zeros<uvec>(0)) ;
  vec    elt_norm  (vec x, uvec pk = zeros<uvec>(0)) ;
  double dual_norm (vec x, uvec pk = zeros<uvec>(0)) ;
  vec proximal(vec x, double lambda, uvec pk = zeros<uvec>(0)) ;

};

#endif
