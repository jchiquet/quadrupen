#ifndef _quadrupen_PENALTY_H
#define _quadrupen_PENALTY_H

#include <RcppArmadillo.h>

#define ZERO 2e-16 // practical zero

using namespace Rcpp;
using namespace arma;

class Penalty {
public: 
  Penalty()        ;
  Penalty(SEXP PK) ;

  void setPenalty(std::string penName) ;

  vec    elt_norm  (vec x) ;
  double pen_norm  (vec x) ;
  double dual_norm (vec x) ;
  vec proximal(vec x, double lambda) ;

private:
  uvec pk ;

  typedef vec (Penalty::*elt_norm_ptr)(vec x) ;
  elt_norm_ptr _current_elt_norm_ptr ;

  typedef double (Penalty::*pen_norm_ptr)(vec x) ;
  pen_norm_ptr _current_pen_norm_ptr ;

  typedef double (Penalty::*dual_norm_ptr)(vec x) ;
  dual_norm_ptr _current_dual_norm_ptr ;

  typedef vec (Penalty::*proximal_ptr)(vec x, double lambda) ;
  proximal_ptr _current_proximal_ptr ;

  vec    elt_norm_L1  (vec x) ;
  double pen_norm_L1  (vec x) ;
  double dual_norm_L1 (vec x) ;
  vec    proximal_L1  (vec x, double lambda) ;

  vec    elt_norm_LINF  (vec x) ;
  double pen_norm_LINF  (vec x) ;
  double dual_norm_LINF (vec x) ;
  vec    proximal_LINF  (vec x, double lambda) ;

  vec    elt_norm_L1L2  (vec x) ;
  double pen_norm_L1L2  (vec x) ;
  double dual_norm_L1L2 (vec x) ;
  vec    proximal_L1L2  (vec x, double lambda) ;

  vec    elt_norm_L1LINF  (vec x) ;
  double pen_norm_L1LINF  (vec x) ;
  double dual_norm_L1LINF (vec x) ;
  vec    proximal_L1LINF  (vec x, double lambda) ;

  vec    elt_norm_COOP  (vec x) ;
  double pen_norm_COOP  (vec x) ;
  double dual_norm_COOP (vec x) ;
  vec    proximal_COOP  (vec x, double lambda) ;
};

#endif
