/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_GROUP_LASSO_MAIN_H
#define _quadrupen_GROUP_LASSO_MAIN__H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "quadrupen_headers.hpp"

using namespace Rcpp;
using namespace arma;

RcppExport SEXP group_lasso(SEXP BETA0    ,
			    SEXP X        , // matrix of features
			    SEXP Y        , // vector of response
			    SEXP TYPE     , // -1="l1/linf", 2="l1/l2" group-Lasso
			    SEXP PK       , // successive groups length
			    SEXP OMEGA    , // structuring matrix
			    SEXP LAMBDA1  ,
			    SEXP NLAMBDA1 ,
			    SEXP MIN_RATIO,
			    SEXP PENSCALE ,
			    SEXP LAMBDA2  ,
			    SEXP INTERCEPT,
			    SEXP NORMALIZE,
			    SEXP WEIGHTS  ,
			    SEXP NAIVE    ,
			    SEXP EPS      ,
			    SEXP MAXITER  ,
			    SEXP MAXFEAT  ,
			    SEXP FUN_OPTIM,
			    SEXP VERBOSE  ,
			    SEXP SPARSE   ) ;

#endif
