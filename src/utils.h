/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_UTILS_H
#define _quadrupen_UTILS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

#define ZERO 2e-16 // practical zero

inline vec get_lambda1(SEXP LAMBDA1, uword n_lambda, double min_ratio, double lmax) {
  vec lambda1 ;
  if (LAMBDA1 != R_NilValue) {
    lambda1  = as<vec>(LAMBDA1)  ;
  } else {
    lambda1 = exp10(linspace(log10(lmax), log10(min_ratio*lmax), n_lambda)) ;
  }
  return(lambda1);
}

vec grp_norm(vec x, uvec pk, int norm, int rep) ;
vec grp_sign(vec x, uvec pk) ;

double get_df_enet(const double &lambda2, mat &R, mat &xAtxA, const sp_mat &S, uvec &A, const uword &fun) ;
double get_df_breg(const double &lambda2, mat &xtx, sp_mat &S, uvec &A) ;

vec  cg(mat A, vec b, vec x, double tol) ;
vec pcg(mat A, mat P, vec b, vec x, double tol) ;

void cholupdate(mat &R, mat& XtX) ;

void choldowndate(mat &R, int j) ;

void bound_to_optimal(vec &betaA, mat &xAtxA, const vec &xty, vec &grd, double &lambda1, const double &lambda2, const double &normy, uvec &A, const uword &monitor, vec &J_hat, vec &D_hat) ;

template <typename any_mat>
void add_var_enet(const uword &n, uword &nbr_in, uword &var_in, vec &betaA, uvec &A, any_mat &x, any_mat &xt, mat &xtxA, mat &xAtxA, mat &xtxw, mat &R, const double &lambda2, const vec &xbar, const sp_mat &spS, const bool &usechol, const uword &fun) {
  
  vec  new_col   ; // column currently added to xtxA
  
  A.resize(nbr_in+1)     ; // update the active set
  A[nbr_in] = var_in     ;
  betaA.resize(nbr_in+1) ; // update the vector of active parameters
  betaA[nbr_in]  = 0.0   ;
  
  new_col = xt * x.col(var_in) - n * xbar * as_scalar(xbar[var_in]);
  if (lambda2 > 0) {
    // Adding the column corresponding to the structurating matrix
    new_col += spS.col(var_in);
  }
  
  // UPDATE THE xtxA AND xAtxA MATRICES
  if (nbr_in > 0) {
    xAtxA = join_cols(xAtxA, xtxA.row(var_in)) ;
  }
  xtxA  = join_rows(xtxA, new_col) ;
  xAtxA = join_rows(xAtxA, trans(xtxA.row(var_in))) ;
  
  if ((fun == 0) & (usechol == 1)) {
    cholupdate(R, xAtxA) ;
  }
  
  if (fun == 1) {
    xtxw.resize(nbr_in+1) ;
    xtxw(nbr_in) = dot(xAtxA.col(nbr_in),betaA);
  }
}


void remove_var_enet(uword &nbr_in, uvec &are_in, vec &betaA, uvec &A, mat &xtxAS, mat &xAtxA, mat &xtxw, mat &R, uvec &null, const bool &usechol, const uword &fun) ;

template <typename any_mat>
void standardize(any_mat &x, const vec &y, const bool &intercept, const bool &normalize, const vec &penscale,
		 vec &xty, vec &normx, double &normy, vec &xbar, double &ybar) {
  
  uword n = x.n_rows;
  uword p = x.n_cols;
  
  if (intercept == 1) {
    xbar = trans(rowvec(mean(x, 0)));
    ybar = mean(y) ;
  } else {
    xbar = zeros(p) ;
    ybar = 0;
  }

  if (normalize == 1) {
    normx = sqrt(trans(sum(square(x),0)) - n * square(xbar));
    for (uword i=0; i<p; i++) {
      x.col(i) /= normx(i);
    }
    xbar /= normx ;
  } else {
    normx = ones(p);
  }
  normy = sqrt(sum(square(y))) ;

  if (any(penscale != 1)) {
    for (uword i=0; i<n; i++) {
       x.row(i) /= trans(penscale) ;
    }
    xbar /= penscale;
  }

  if (intercept == 1) {
    xty = trans(trans(y-ybar) * x) ;
    for (uword i=0; i<p; i++) {
       xty(i) -=  sum(y-ybar) * xbar(i);
    }
  } else {
    xty = trans(y.t()*x) ;
  }
}

#endif

