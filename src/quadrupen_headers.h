/***
 * @file interface.h
 *
 */
#ifndef __QUADRUPEN_HEADERS_HPP
#define __QUADRUPEN_HEADERS_HPP

// Include Armadillo / Rcpp / R to C/C++ basics
#include <string.h>
#include <sys/time.h>
#include <RcppArmadillo.h>

#define ZERO 2e-16 // practical zero

// Include utils and optimization routines
#include "utils.h"
#include "quadratic.h"
#include "first_order.h"
// #include "data_reg.h"
// #include "path.h"
// #include "active_set.h"
// #include "penalties.h"

#endif

