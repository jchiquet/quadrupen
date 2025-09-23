/***
 * @file interface.hpp
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
// #include "data_reg.hpp"
// #include "path.hpp"
// #include "active_set.hpp"
// #include "penalties.hpp"

#endif

