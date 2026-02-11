/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#ifndef _quadrupen_QUADRATIC_H
#define _quadrupen_QUADRATIC_H

#define ARMA_WARN_LEVEL 1

#include <RcppArmadillo.h>
#include "utils.h"

using namespace Rcpp;
using namespace arma;

int quadra_enet(vec& x0, mat& R,  mat& xAtxA, vec xty, vec sgn_grd, double &pen, uvec& null, bool usechol, double tol) ;

#endif
