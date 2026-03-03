/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/GroupElasticNet.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List group_enet_l1l2_dense_cpp(
    const Environment &dataModel, // data structure
    const bool        &intercept, // Boolean for intercept
    const arma::uvec  &group    , // Vector of groups
    const List        &regParam , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  GroupElasticNet<mat> grpenet(data, group, regParam, control);
  
  List results = grpenet.solution_path(control);

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = grpenet.path_tuning(),
      Named("l2") = grpenet.gamma_
    ),
    Named("coef")               = grpenet.coefficients(),
    Named("coef_debiased")      = grpenet.debiased_coefficients(),
    Named("active")             = grpenet.active_var(),
    Named("intercept")          = grpenet.intercept_,
    Named("intercept_debiased") = grpenet.intercept_debiased_,
    Named("normx")              = grpenet.data_.norm_X_,
    Named("df")                 = grpenet.df_,
    Named("monitoring")         = results
  );
  
}

// [[Rcpp::export]]
List group_enet_l1l2_sparse_cpp(
    const Environment &dataModel, // data structure
    const bool        &intercept, // Boolean for intercept
    const arma::uvec  &group    , // Vector of groups
    const List        &regParam , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  
  GroupElasticNet<sp_mat> grpenet(data, group, regParam, control);
  
  List results = grpenet.solution_path(control);
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = grpenet.path_tuning(),
      Named("l2") = grpenet.gamma_
    ),
    Named("coef")               = grpenet.coefficients(),
    Named("coef_debiased")      = grpenet.debiased_coefficients(),
    Named("active")             = grpenet.active_var(),
    Named("intercept")          = grpenet.intercept_,
    Named("intercept_debiased") = grpenet.intercept_debiased_,
    Named("normx")              = grpenet.data_.norm_X_,
    Named("df")                 = grpenet.df_,
    Named("monitoring")         = results
  );
  
}
