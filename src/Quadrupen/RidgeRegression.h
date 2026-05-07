/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_RIDGE_H
#define _quadrupen_RIDGE_H

#include "Regularizer.h"
#include "PenaltyDense.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

class RidgeRegression : public Regularizer<mat>{
public:
  RidgeRegression(RegressionData<mat>&, const List&);

  double get_lambda_max() {return(penalty_.dual_norm(data_.XTy_, lambda_factor_));}
  
  List solution_path(const mat&);
  
  DensePenalty<DenseNorm::L2> penalty_ ; // main penalty object 
  
};

#endif

