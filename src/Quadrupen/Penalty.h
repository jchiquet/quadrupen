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

enum class SimpleNorm {L1, LINF, RIDGE};
enum class MixedNorm {L1L2, L1LINF, COOP};

using namespace Rcpp;
using namespace arma;
using namespace std;

template <SimpleNorm norm> class SimplePenalty {
public: 
  
  SimplePenalty() {} ;

  double pen_norm  (vec x) ;
  vec    elt_norm  (vec x) ;
  double dual_norm (vec x) ;
  vec proximal(vec x, double lambda) ;

};

template <MixedNorm norm> class MixedPenalty {
public: 
  
  MixedPenalty() {} ;
  
  vec    elt_norm  (vec x, uvec pk) ;
  double pen_norm  (vec x, uvec pk, vec wk) ;
  double dual_norm (vec x, uvec pk, vec wk) ;
  vec proximal(vec x, vec lambda, uvec pk) ;
  
};

#endif
