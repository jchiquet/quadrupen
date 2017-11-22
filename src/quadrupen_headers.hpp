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
#include "utils/utils.h"
#include "optimization/quadratic.h"
#include "optimization/first_order.h"
#include "tmp_cpp/data_reg.hpp"
#include "tmp_cpp/path.hpp"
#include "tmp_cpp/active_set.hpp"
#include "tmp_cpp/penalties.hpp"

#endif

