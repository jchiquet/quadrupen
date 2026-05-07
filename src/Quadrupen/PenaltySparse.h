#ifndef _quadrupen_PENALTY_SPARSE_H
#define _quadrupen_PENALTY_SPARSE_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>

enum class SparseNorm {L1, MCP, SCAD};

using namespace Rcpp;
using namespace arma;
using namespace std;

template <SparseNorm norm> class SparsePenalty {
  public: 
    
    double gamma_ = 0 ; // Optional factor for non-concave penalties
    SparsePenalty() {} ;
    
    vec    elt_norm  (const vec& x, const vec& w, double lambda=0) ;
    vec    elt_dual_norm  (const vec& x, const vec& w, double lambda=0) ;
    double pen_norm  (const vec& x, const vec& w, double lambda=0) ;
    double dual_norm (const vec& x, const vec& w, double lambda=0) ;
    vec proximal(const vec& x, double lambda, const vec& w) ;
    double lambda_max (const vec& XTy, const vec& w) ;
    vec optimality(const vec& x, double lambda, const vec& w)  ;
    
};

template<SparseNorm norm>
  vec SparsePenalty<norm>::optimality(const vec& grad, double lambda, const vec& w)  {
    return(elt_dual_norm(grad, w) - lambda) ;
  }

template<SparseNorm norm>
  double SparsePenalty<norm>::lambda_max(const vec& XTy, const vec& w)  {
    return(dual_norm(XTy, w)) ;
  }

#endif
