/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _BoundedRegression_H
#define _BoundedRegression_H

#include "Regularizer.h"
#include "ActiveSet.h"
#include "Penalty.h"
#include "OptimizerLINF.h"

#define ZERO 2e-16 // practical zero

using namespace Rcpp;
using namespace arma;
using namespace std;

class BoundedRegression : public Regularizer<mat> {
public:
  
  BoundedRegression(RegressionData<mat>&, const List&, const List&);

  double get_lambda_max() {return(penalty_.dual_norm(data_.XTy_));}
  
  List solution_path(const List&);

  // Specific to Bounded regression
  SimplePenalty<SimpleNorm::LINF> penalty_ ; // main penalty object 
  OptimizerLINF<mat> solver_ ; // Solvers for LINF penalty
  ActiveSet<mat> set_       ; // Active set of variable and data

  // Compute degrees of freedom for the current estimate
  double get_df() ; 
  
};

#endif

