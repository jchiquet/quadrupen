#ifndef _quadrupen_PENALTY_SIMPLE_H
#define _quadrupen_PENALTY_SIMPLE_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>

#define ZERO 2e-16 // practical zero

enum class SimpleNorm {L1, LINF, RIDGE};

using namespace Rcpp;
using namespace arma;
using namespace std;

template <SimpleNorm norm> class SimplePenalty {
public: 
  
  SimplePenalty() {} ;

  vec    elt_norm  (const vec& x) ;
  vec    elt_dual_norm  (const vec& x) ;
  double pen_norm  (const vec& x, const vec& w) ;
  double dual_norm (const vec& x, const vec& w) ;
  vec proximal(const vec& x, double lambda, const vec& w) ;

};

#endif
