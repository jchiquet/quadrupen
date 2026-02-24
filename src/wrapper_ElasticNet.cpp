/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/ElasticNet.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List elastic_net_dense_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
    ) {

  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  ElasticNet<mat> enet(data, regParam, control);

  List results = enet.solution_path(control);

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = enet.path_tuning(),
      Named("l2") = enet.gamma_
    ),
    Named("coef")               = enet.coefficients(),
    Named("coef_debiased")      = enet.debiased_coefficients(),
    Named("active")             = enet.active_var(),
    Named("intercept")          = enet.intercept_,
    Named("intercept_debiased") = enet.intercept_debiased_,
    Named("normx")              = enet.data_.norm_X_,
    Named("df")                 = enet.df_,
    Named("monitoring")         = results
  );
  
}

// [[Rcpp::export]]
List elastic_net_sparse_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  
  ElasticNet<sp_mat> enet(data, regParam, control);
  
  List results = enet.solution_path(control);
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = enet.path_tuning(),
      Named("l2") = enet.gamma_
    ),
    Named("coef")               = enet.coefficients(),
    Named("coef_debiased")      = enet.debiased_coefficients(),
    Named("active")             = enet.active_var(),
    Named("intercept")          = enet.intercept_,
    Named("intercept_debiased") = enet.intercept_debiased_,
    Named("normx")              = enet.data_.norm_X_,
    Named("df")                 = enet.df_,
    Named("monitoring")         = results
  );
  
}
