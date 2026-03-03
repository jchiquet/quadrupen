/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_GENERIC_OPTIMIZER_H
#define _quadrupen_GENERIC_OPTIMIZER_H

#include "RegressionData.h"
#include "ActiveSet.h"
#include "ActiveSetGroup.h"
#include "Penalty.h"

enum SolverType {FISTA, QUADRA};

#define ZERO 2e-16 // practical zero

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix> class Optimizer {
public:
  
  Optimizer() {} ;
  Optimizer(const List& control) ;
  
  SolverType algorithm_ ;
  bool verbosity_       ;
  uword iter_, maxiter_, maxfeat_, monitoring_ ;
  double accuracy_, gap_, J_, D_ ; 
  vector<uword> inner_iter_   ; 
  vector<double> J_vec_, D_vec_ ;   
  
  uword conjugate_gradient(
      vec& x0,
      const mat& A,
      const vec& b,
      const double& accuracy,
      const uword& max_iter) ;

};

template <typename matrix>
Optimizer<matrix>::Optimizer(const List& control) :
  accuracy_(control["threshold"]),
  verbosity_(control["verbose"]),
  maxiter_(control["maxiter"]),
  maxfeat_(control["maxfeat"]),
  monitoring_(control["monitor"]) {
  
  if (as<std::string>(control["method"]) == "FISTA") algorithm_ = FISTA;
  if (as<std::string>(control["method"]) == "QUADRA") algorithm_ = QUADRA;
  
}

template <typename matrix>
uword Optimizer<matrix>::conjugate_gradient(
    vec& x0,
    const mat& A,
    const vec& b,
    const double& accuracy,
    const uword& max_iter) {
  
  vec r = b - A * x0;
  vec p = r ;
  double rs_old = sum(square(r)) ;
  
  double rs_new = rs_old ;
  uword i = 0;
  double alpha ;
  mat Ap ;
  
  while ((sqrt(rs_new) > accuracy) & (i < max_iter)) {
    Ap = A * p;
    alpha = rs_old/dot(p,Ap) ;
    if (std::isfinite(alpha)) {
      x0 += alpha * p ;
      r -= alpha * Ap ;
      // Polak-Ribière update
      rs_new = dot(r,-alpha * Ap);
      p = r + rs_new/rs_old*p;
      rs_old = rs_new;
      i++;
    } else {
      break ;
    }
  }
  return(i) ;
}

#endif

