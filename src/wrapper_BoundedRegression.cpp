/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/BoundedRegression.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List bounded_regression_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
    ) {

  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  BoundedRegression bounded(data, regParam, control);

  List results = bounded.solution_path(control);

  return List::create(
    Named("tuning_param") = List::create(
      Named("linf") = bounded.lambdas_,
      Named("l2")   = bounded.gamma_
    ),
    Named("coef")        = bounded.coefficients(),
    Named("active")      = bounded.unbounded_var(),
    Named("intercept")   = bounded.intercept_,
    Named("normx")       = bounded.data_.norm_X_,
    Named("df")          = bounded.df_,
    Named("monitoring")  = results
  );

}
