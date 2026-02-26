/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_RIDGE_H
#define _quadrupen_RIDGE_H

#include "Regularizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

class RidgeRegression : public SimpleRegularizer<mat,SimpleNorm::RIDGE>{
public:
  RidgeRegression(RegressionData<mat>&, const List&);
  
  double get_lambda_max() {return(penalty_.dual_norm(data_.XTy_));}
  
  List solution_path(const mat&);
  
};

#endif

