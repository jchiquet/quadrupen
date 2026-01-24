/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "RidgeRegression_class.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List ridge_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const List        &regParam    , // regularization parameters
    const List        &controlFit    // config of the optimization 
    ) {

  RegressionData<mat> data(dataModel, intercept, as<bool>(controlFit["normalize"])) ;

  RidgeRegression ridge(data, intercept, regParam);

  List results = ridge.solution_path(trimatu(as<mat>(dataModel["C_inv"])));
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l2") = ridge.tuning_param(),
      Named("l1") = 0 
    ),
    Named("beta")        = ridge.coefficients(),
    Named("mu")          = ridge.interceptTerm(),
    Named("normx")       = ridge.data().norm_X(),
    Named("df")          = ridge.degrees_freedom(),
    Named("monitoring")  = results
  );
  
}
