/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "ElasticNet_class.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List elastic_net_dense_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
    ) {

  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  ElasticNet enet(data, intercept, regParam, control);

  List results = enet.solution_path(control);
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = enet.path_tuning(),
      Named("l2") = enet.struct_tuning()
    ),
    Named("beta")        = enet.coefficients(),
    Named("active")      = enet.active_var(),
    Named("mu")          = enet.interceptTerm(),
    Named("normx")       = enet.data().norm_X(),
    Named("df")          = enet.degrees_freedom(),
    Named("monitoring")  = results
  );
  
}
