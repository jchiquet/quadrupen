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

  void optimality_gap(
      vec& beta,
      vec& grad,
      const double& lambda,
      const double& gamma,
      RegressionData<matrix> &data,
      ActiveSet<matrix>& set, uword type) ;
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

template <typename matrix> 
void Optimizer<matrix>::optimality_gap(
    vec& beta,
    vec& grad,
    const double& lambda,
    const double& gamma,
    RegressionData<matrix> &data,
    ActiveSet<matrix>& set, uword type) {
  
  // nu equals the max |gradient|
  double nu = arma::norm(grad, "inf");
  double loss = .5 * pow(data.norm_y_ ,2) + 
    dot(beta, .5 * set.XATXA_ * beta - data.XTy_(set.A_)) ;
  double old_J = J_, old_D = D_ ;
  J_ = loss - dot(beta, grad(set.A_))  ;
  uvec Ac ;
  
  switch (type) {
  case 1: // Grandvalet's bound
    Ac = find(grad > nu); // set of adversarial variables outside the boundary
    D_ = J_ * (1 - lambda/nu) - 
      (pow(lambda,2)/(2*gamma))*((lambda*(data.p_-Ac.n_elem))/nu + 
      pow(arma::norm(grad(Ac),2)/nu,2)-data.p_);
    break;
  case 2: // Fenchel's bound
    if (nu < lambda) nu = lambda;
    D_ = loss * (1+pow(lambda/nu,2)) + sum(abs(lambda*beta)) + 
      (lambda/nu)*(dot(beta,data.XTy_(set.A_))-pow(data.norm_y_,2));
    break;
  default: 
    D_ = datum::inf ;
  break;
  }
  
  // keep the smallest bound reached so far for a given lambda value
  if ((old_J < J_) && (old_D - D_) < (old_J - J_)) {D_ = old_D ; }
  
}

#endif

