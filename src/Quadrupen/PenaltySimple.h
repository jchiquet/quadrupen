#ifndef _quadrupen_PENALTY_SIMPLE_H
#define _quadrupen_PENALTY_SIMPLE_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>

enum class SimpleNorm {L1, L2, LINF, MCP, SCAD};

using namespace Rcpp;
using namespace arma;
using namespace std;

template <SimpleNorm norm> class SimplePenalty {
public: 
  
  double gamma_ = 0 ; // Optional factor for non-concave penalties
  SimplePenalty() {} ;

  vec    elt_norm  (const vec& x, const vec& w, double lambda=0) ;
  vec    elt_dual_norm  (const vec& x, const vec& w, double lambda=0) ;
  double pen_norm  (const vec& x, const vec& w, double lambda=0) ;
  double dual_norm (const vec& x, const vec& w, double lambda=0) ;
  vec proximal(const vec& x, double lambda, const vec& w) ;
  double lambda_max (const vec& XTy, const vec& w) ;
  vec optimality(const vec& x, double lambda, const vec& w)  ;
  
};

template<SimpleNorm norm>
vec SimplePenalty<norm>::optimality(const vec& grad, double lambda, const vec& w)  {
  return(elt_dual_norm(grad, w) - lambda) ;
}

template<SimpleNorm norm>
double SimplePenalty<norm>::lambda_max(const vec& XTy, const vec& w)  {
  return(dual_norm(XTy, w)) ;
}

#endif
