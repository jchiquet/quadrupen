#ifndef _quadrupen_PENALTY_GROUP_H
#define _quadrupen_PENALTY_GROUP_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>

#define ZERO 2e-16 // practical zero

enum class GroupNorm {L1L2, L1LINF, COOP};

using namespace Rcpp;
using namespace arma;
using namespace std;

template <GroupNorm norm> class GroupPenalty {
public: 
  
  double alpha_ = 0 ; // Optional sparse factor 
  GroupPenalty() {} ;
  GroupPenalty(double alpha) : alpha_(alpha){} ;
  
  vec    elt_norm  (const vec& x, const uvec& pk, const vec& wk)         ;
  vec    elt_dual_norm  (const vec& x, const uvec& pk, const vec& wk)    ;
  double pen_norm  (const vec& x, const uvec& pk, const vec& wk) ;
  double dual_norm (const vec& x, const uvec& pk, const vec& wk) ;
  vec proximal(const vec& x, double lambda, const uvec& pk, const vec& wk)  ;
  
};

template<GroupNorm norm>
double GroupPenalty<norm>::pen_norm(const vec& x, const uvec& pk, const vec& wk) {
  return(sum(elt_norm(x, pk, wk)));
}

template<GroupNorm norm>
double GroupPenalty<norm>::dual_norm(const vec& x, const uvec& pk, const vec& wk) {
  return(max(elt_dual_norm(x, pk, wk)));
}

#endif
