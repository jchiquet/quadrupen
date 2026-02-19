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
  
  // Data declaration and standardization
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
  vec lambda_factor = as<vec>(regParam["lambda_factor"]) ;
  double gamma = as<double>(regParam["gamma"]) ;
  data.scale_struct(sqrt(data.n_*gamma)*pow(lambda_factor,-1)) ;
  data.scale_regressors(lambda_factor) ;
  
  // Create the scaled/transformed data
  mat C_inv = solve(trimatu(chol(data.S_.as_dense())), eye(data.p_, data.p_)) ;
  mat U, V ; vec D ; // U D V = X C^-1
  svd_econ(U, D, V, (data.X_.each_row() - data.X_bar_.t())*C_inv) ;
  mat Proj = U * diagmat(square(D)/(square(D) + 1)) * U.t()  ;
  mat K12 = sqrtmat_sympd(diagmat(ones(data.n_)) - Proj) ;
  mat X_tilde = K12 * (data.X_.each_row() - data.X_bar_.t()) ;
  vec y_tilde = K12 * (data.y_ - data.y_bar_) ;
  
  // Create the corresponding scaled/transformed data 
  RegressionData<mat> scaled_data(
      X_tilde,
      y_tilde,
      sp_mat(data.p_, data.p_),
      ones(data.p_), false, false) ;

  Lava<mat> lava(scaled_data, Proj, regParam, control);

  List results = lava.solution_path(control);

  lava.post_treatment(data, U, V, D, C_inv) ;
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = lava.path_tuning(),
      Named("l2") = gamma
    ),
    Named("beta")        = lava.sparse_coef_,
    Named("b")           = lava.dense_coef_,
    Named("active")      = lava.active_var(),
    Named("mu")          = lava.intercept(),
    Named("normx")       = lava.data().norm_X_,
    Named("df")          = lava.degrees_freedom(),
    Named("monitoring")  = results
  );
  
}

// // [[Rcpp::export]]
// List lava_sparse_cpp(
//     const Environment &dataModel   , // data structure
//     const bool        &intercept   , // Boolean for intercept
//     const List        &regParam    , // regularization parameters
//     const List        &control    // config of the optimization 
// ) {
//   
//   RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
//   
//   Lava<sp_mat> lava(data, intercept, regParam, control);
//   
//   List results = lava.solution_path(control);
//   
//   return List::create(
//     Named("tuning_param") = List::create(
//       Named("l1") = lava.path_tuning(),
//       Named("l2") = lava.struct_tuning()
//     ),
//     Named("beta")        = lava.coefficients(),
//     Named("active")      = lava.active_var(),
//     Named("mu")          = lava.intercept(),
//     Named("normx")       = lava.data().norm_X_,
//     Named("df")          = lava.degrees_freedom(),
//     Named("monitoring")  = results
//   );
//   
// }
