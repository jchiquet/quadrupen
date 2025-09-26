/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_PATH_H
#define _quadrupen_PATH_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

class PATH {

protected:
  vec penLevels    ; // the successive levels of penalties along the path
  uword maxFeat    ; // maximal number of features to enter the path

  // VARIABLES INDEXED BY PENALTY LEVELS
  sp_mat beta   ; // final matrix of fitted coefficients
  vec mu        ; // vector of the successive intercept term
  vec df        ; // vector of degrees of freedom along the path
  mat  nonzeros ; // contains non-zero value of beta
  mat  iA       ; // contains row indices of the non-zero values
  mat  jA       ; // contains column indices of the non-zero values

public:
  PATH(SEXP MAXFEAT) ;

  void grid_penLevels(SEXP PENLEVELS,
		      SEXP PENLEN   ,
		      SEXP MINRATIO ,
		      double LMAX ) ;

  const vec   & get_penLevels() const { return penLevels ; }
  const uword & get_maxFeat()   const { return maxFeat   ; }
  const vec   & get_mu()        const { return mu ; }
  const vec   & get_df()        const { return df ; }

};

#endif
