/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _ElasticNet_H
#define _ElasticNet_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "GenericRegularizer_class.h"
#include "ActiveSet_class.h"
#include "Optimizer_class.h"

#define ZERO 2e-16 // practical zero

using namespace Rcpp;
using namespace arma;
using namespace std;

class ElasticNet : public GenericRegularizer<mat> {

  public:
  ElasticNet(RegressionData<mat>&, const bool&, const List&, const List&);
  
  List solution_path(const List&);

private:

  // Specific to Elastic-Net regularization
  double gamma_  ; // overall amount of l2 penalty

  mat  R_     ; // Cholesky decomposition of XAtXA
  vec  xtxw_  ; // t(x_A) * x_A * beta(A)
  urowvec iA_ ; // contains row indices of the non-zero values
  urowvec jA_ ; // contains column indices of the non-zero values

  // Compute degrees of freedom for the current estimate
  double get_df() ; 

};

#endif

