/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/Lava.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List lava_dense_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  
  Lava<mat> lava(data, intercept, regParam, control);
  
  List results = lava.solution_path(control);
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = lava.path_tuning(),
      Named("l2") = lava.struct_tuning()
    ),
    Named("delta")       = lava.sparse_coefficients(),
    Named("beta")        = lava.dense_coefficients(),
    Named("active")      = lava.active_var(),
    Named("mu")          = lava.interceptTerm(),
    Named("normx")       = lava.data().norm_X_,
    Named("df")          = lava.degrees_freedom(),
    Named("monitoring")  = results
  );
  
}

// [[Rcpp::export]]
List lava_sparse_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  
  Lava<sp_mat> lava(data, intercept, regParam, control);
  
  List results = lava.solution_path(control);
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = lava.path_tuning(),
      Named("l2") = lava.struct_tuning()
    ),
    Named("beta")        = lava.coefficients(),
    Named("active")      = lava.active_var(),
    Named("mu")          = lava.interceptTerm(),
    Named("normx")       = lava.data().norm_X_,
    Named("df")          = lava.degrees_freedom(),
    Named("monitoring")  = results
  );
  
}
