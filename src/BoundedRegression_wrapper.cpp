/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "BoundedRegression_class.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List bounded_regression_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
    ) {

  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  BoundedRegression bounded(data, intercept, regParam, control);

  List results = bounded.solution_path(control);

  return List::create(
    Named("tuning_param") = List::create(
      Named("linf") = bounded.path_tuning(),
      Named("l2")   = bounded.struct_tuning()
    ),
    Named("beta")        = bounded.coefficients(),
    Named("mu")          = bounded.interceptTerm(),
    Named("normx")       = bounded.data().norm_X(),
    Named("df")          = bounded.degrees_freedom(),
    Named("monitoring")  = results
  );

}
