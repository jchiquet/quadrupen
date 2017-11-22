/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_FIRSTORDER_H
#define _quadrupen_FIRSTORDER_H

#include <RcppArmadillo.h>
#include "utils/utils.h"

using namespace Rcpp;
using namespace arma;

#define ZERO 2e-16 // practical zero

int fista_lasso(vec &x0, mat &xtx, vec xty, double &pen, uvec &null, double &L0, double eps) ;
int fista_breg (vec &x0, const mat &xtx, const vec &xty, vec& grd, double &pen, double &L0, double eps) ;
int fista_grp(vec &x0, uvec &pk, mat &xtx, vec xty, vec &pen, uword &grpnorm, double &L0, double eps) ;

vec proximal_grp2(vec u, uvec pk, vec lambda, double L) ;
vec proximal_grpinf(vec u, uvec pk, vec lambda, double L) ;
vec proximal_inf(vec v, vec lambda) ;

int pathwise_enet(vec& x0, mat& xtx, vec xty, vec& xtxw, double& pen, uvec &null, double& gam, double eps) ;

// int pathwise_grpl1linf(vec &x0, uvec &pk, mat &xtx, vec xty, vec &xtxw, vec &pen, double &gam, double eps) ;


#endif
