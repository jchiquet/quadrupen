/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_BOUNDEDREG_H
#define _quadrupen_BOUNDEDREG_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "quadrupen_headers.hpp"

RcppExport SEXP bounded_reg(SEXP X        ,
			    SEXP Y        ,
			    SEXP STRUCT   ,
			    SEXP LAMBDA1  ,
			    SEXP N_LAMBDA ,
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
			    SEXP FUN      ,
			    SEXP VERBOSE  ,
			    SEXP SPARSE   ,
			    SEXP BULLETPROOF) ;

#endif
