/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _BoundedRegression_H
#define _BoundedRegression_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "GenericRegularizer_class.h"
#include "Optimizer_class.h"

#define ZERO 2e-16 // practical zero

using namespace Rcpp;
using namespace arma;
using namespace std;

class BoundedRegression : public GenericRegularizer {
public:
  BoundedRegression(RegressionData&, const bool&, const List&);
  
  List solution_path(const List&);

private:

  // Specific to Bounded regression
  mat    XTX        ; // Gram matrix
  double gamma_     ; // overall amount of l2 penalty
  sp_mat S_scaled   ; // structuring matrix scaled according to main penalty factors
  uvec B            ; // guys reaching the boundary
  urowvec iB_       ; // contains row indices of the bounded variables
  urowvec jB_       ; // contains column indices of the bounded variables
  
  // Compute degrees of freedom for the current estimate
  double get_df() ; 
  
};

#endif

