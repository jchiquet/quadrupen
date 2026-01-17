/***
 * @file interface.h
 *
 */
#ifndef __QUADRUPEN_HEADERS_HPP
#define __QUADRUPEN_HEADERS_HPP

#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR

// [[Rcpp::depends(RcppArmadillo)]]

// Include Armadillo / Rcpp / R to C/C++ basics
#include "RcppArmadillo.h"
#include <string.h>
#include <sys/time.h>

#define ZERO 2e-16 // practical zero

// Include utils and optimization routines
#include "utils.h"
#include "quadratic.h"
#include "first_order.h"

#endif

