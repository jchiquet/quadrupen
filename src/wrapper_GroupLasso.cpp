/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/GroupLasso.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List group_lasso_l1l2_dense_cpp(
    const Environment &dataModel, // data structure
    const bool        &intercept, // Boolean for intercept
    const arma::uvec  &group    , // Vector of groups
    const List        &regParam , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  GroupLassoL1L2<mat> grplasso(data, group, regParam, control);
  
  List results = grplasso.solution_path(control);

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = grplasso.path_tuning(),
      Named("l2") = grplasso.gamma_
    ),
    Named("coef")               = grplasso.coefficients(),
    Named("coef_debiased")      = grplasso.debiased_coefficients(),
    Named("active")             = grplasso.active_var(),
    Named("intercept")          = grplasso.intercept_,
    Named("intercept_debiased") = grplasso.intercept_debiased_,
    Named("normx")              = grplasso.data_.norm_X_,
    Named("df")                 = grplasso.df_,
    Named("monitoring")         = results
  );
  
}

// [[Rcpp::export]]
List group_lasso_l1l2_sparse_cpp(
    const Environment &dataModel, // data structure
    const bool        &intercept, // Boolean for intercept
    const arma::uvec  &group    , // Vector of groups
    const List        &regParam , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  
  GroupLassoL1L2<sp_mat> grplasso(data, group, regParam, control);
  
  List results = grplasso.solution_path(control);
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = grplasso.path_tuning(),
      Named("l2") = grplasso.gamma_
    ),
    Named("coef")               = grplasso.coefficients(),
    Named("coef_debiased")      = grplasso.debiased_coefficients(),
    Named("active")             = grplasso.active_var(),
    Named("intercept")          = grplasso.intercept_,
    Named("intercept_debiased") = grplasso.intercept_debiased_,
    Named("normx")              = grplasso.data_.norm_X_,
    Named("df")                 = grplasso.df_,
    Named("monitoring")         = results
  );
  
}
