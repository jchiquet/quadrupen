/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_GENERIC_OPTIMIZER_H
#define _quadrupen_GENERIC_OPTIMIZER_H

#include "RegressionData.h"
#include <functional>

enum SolverType {FISTA, QUADRA, PGD};

using namespace Rcpp;
using namespace arma;
using namespace std;

class Optimizer {
public:
  
  Optimizer() {} ;
  Optimizer(const List& control) ;

  SolverType algorithm_ ;
  double accuracy_, gap_, J_, D_ ; 
  bool verbosity_       ;
  uword iter_, maxiter_, maxfeat_, monitoring_ ;
  vector<uword> inner_iter_   ; 
  vector<double> J_vec_, D_vec_ ;   
  
  uword conjugate_gradient(
      vec& x0,
      const mat& A,
      const vec& b,
      const double& accuracy,
      const uword& max_iter) ;

  double estimate_lipschitz(
      const mat& XTX,
      uword max_it = 15,
      double tol = 1e-4
  ) ;
  
  uword fista(
      vec& beta,
      const double& lambda, 
      const vec& XTy,
      const mat& XTX,
      std::function<vec(const vec&, double)> proximal_operator, 
      const double& accuracy, 
      const uword& max_iter
  ) ;

  uword pgd(
      vec& beta,
      const double& lambda,
      const vec& XTy,
      const mat& XTX,
      std::function<vec(const vec&, double)> proximal_operator, 
      const double& accuracy,
      const uword& max_iter,
      const uword m = 3
  ) ;

  void optimality_gap(
      vec& beta,
      vec& grad,
      const double& lambda,
      const double& gamma,
      const vec& XTy,
      const mat& XTX,
      const double& norm_y,
      uvec A,
      uword type
    ) ;
};

#endif
