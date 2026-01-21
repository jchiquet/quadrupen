/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_OPTIMIZER_H
#define _quadrupen_OPTIMIZER_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "RegressionData_class.h"
#include "Penalty_class.h"
#include "utils.h"

#define ZERO 2e-16 // practical zero

enum SolverType {FISTA, QUADRA};

using namespace Rcpp;
using namespace arma;
using namespace std;

class Optimizer {
public:

  Optimizer() ;
  Optimizer(RegressionData&, Penalty&, SolverType&) ;
  
  int run(vec&, uvec&, vec&, mat&, double&, const double&, const uword&) ;
  
private:
  RegressionData data_ ;
  Penalty penalty_ ;
  SolverType solver_ ;

  typedef int (Optimizer::*solver_ptr)(
      vec& beta,
      uvec& B,
      vec& grad,
      mat& XTX,
      double& lambda,
      const double& accuracy,
      const uword& max_iter) ;
  
  solver_ptr _current_solver_ptr ;

  int quadratic_breg(
      vec& beta,
      uvec& B,
      vec& grad,
      mat& XTX,
      double& lambda,
      const double& accuracy=1e-10,
      const uword& max_iter=50) ;
  
  int fista_breg(
      vec& beta0,
      uvec& B,
      vec& grad,
      mat& XTX, 
      double& lambda, 
      const double& eps=1e-2, 
      const uword& max_iter=10000) ;
  
};

#endif

