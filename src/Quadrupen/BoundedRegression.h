/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _BoundedRegression_H
#define _BoundedRegression_H

#include "GenericRegularizer.h"
#include "OptimizerLINF.h"

#define ZERO 2e-16 // practical zero

using namespace Rcpp;
using namespace arma;
using namespace std;

class BoundedRegression : public GenericRegularizer<mat,Norm::LINF> {
public:
  
  BoundedRegression(RegressionData<mat>&, const List&, const List&);
  
  List solution_path(const List&);

  // Specific to Bounded regression
  OptimizerLINF<mat> solver_ ; // Solvers for LINF penalty
  mat    XTX    ; // Gram matrix
  double gamma_ ; // overall amount of l2 penalty
  vec beta_     ; // vector of current parameters
  vec grad_     ; // vector of current gradient (smooth part)
  urowvec iA_   ; // contains row indices of the active variables
  urowvec jA_   ; // contains column indices of the active variables
  
  // Compute degrees of freedom for the current estimate
  double get_df() ; 
  
};

#endif

