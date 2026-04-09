/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/GroupLava.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List group_lava_l1l2_dense_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const arma::uvec  &group    , // Vector of groups
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  // Data declaration and standardization
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
  vec lambda_factor = as<vec>(regParam["lambda_factor"]) ;
  double gamma = as<double>(regParam["gamma"]) ;
  data.scale_struct(sqrt(gamma)*ones(data.p_)) ;
  // data.scale_struct(sqrt(gamma)*pow(lambda_factor,-1)) ;
  // data.scale_regressors(lambda_factor) ;
  
  // Create the scaled/transformed data
  mat C_inv = solve(trimatu(chol(data.S_.as_dense())), eye(data.p_, data.p_)) ;
  mat U, V ; vec D ; // U D V = X C^-1
  svd_econ(U, D, V, (data.X_.each_row() - data.X_bar_.t())*C_inv) ;
  mat Proj = U * diagmat(square(D)/(square(D) + 1)) * U.t() ;
  mat K12 = U * diagmat(1/sqrt(square(D) + 1)) * U.t() ;
  mat X_tilde = K12 * (data.X_.each_row() - data.X_bar_.t()) ;
  vec y_tilde = K12 * (data.y_ - data.y_bar_) ;
  
  // Create the corresponding scaled/transformed data 
  RegressionData<mat> scaled_data(
      X_tilde,
      y_tilde,
      sp_mat(data.p_, data.p_),
      ones(data.p_), false, false) ;
  
  GroupLava<mat,MixedNorm::L1L2> group_lava(scaled_data, Proj, group, regParam, control);

  List results = group_lava.solution_path(control);

  // Un-normalized vector of dense coefficients
  mat b = (C_inv * V * diagmat(D/(square(D) + 1)) * U.t()) * ( 
    (data.y_ - data.y_bar_) * ones(1, group_lava.lambdas_.size()) - 
      (data.X_.each_row() - data.X_bar_.t()) * group_lava.coefficients()
  ) ;
  
  group_lava.post_treatment(data, b) ;

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = group_lava.lambdas_,
      Named("l2") = gamma
    ),
    
    Named("coef")          = group_lava.coef_,
    Named("coef_debiased") = diagmat(1/data.norm_X_) * b + group_lava.debiased_coefficients(),
    Named("sparse_coef")   = group_lava.sparse_coef_,
    Named("sparse_coef_debiased") = group_lava.debiased_coefficients(),
    Named("active")        = group_lava.active_var(),
    Named("intercept")     = group_lava.intercept_,
    Named("intercept_debiased") = group_lava.intercept_debiased_,
    Named("normx")         = group_lava.data_.norm_X_,
    Named("df")            = group_lava.df_,
    Named("monitoring")    = results
  
  );
  
}

// [[Rcpp::export]]
List group_lava_l1linf_dense_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const arma::uvec  &group    , // Vector of groups
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  // Data declaration and standardization
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
  vec lambda_factor = as<vec>(regParam["lambda_factor"]) ;
  double gamma = as<double>(regParam["gamma"]) ;
  data.scale_struct(sqrt(gamma)*ones(data.p_)) ;
  // data.scale_struct(sqrt(gamma)*pow(lambda_factor,-1)) ;
  // data.scale_regressors(lambda_factor) ;
  
  // Create the scaled/transformed data
  mat C_inv = solve(trimatu(chol(data.S_.as_dense())), eye(data.p_, data.p_)) ;
  mat U, V ; vec D ; // U D V = X C^-1
  svd_econ(U, D, V, (data.X_.each_row() - data.X_bar_.t())*C_inv) ;
  mat Proj = U * diagmat(square(D)/(square(D) + 1)) * U.t() ;
  mat K12 = U * diagmat(1/sqrt(square(D) + 1)) * U.t() ;
  mat X_tilde = K12 * (data.X_.each_row() - data.X_bar_.t()) ;
  vec y_tilde = K12 * (data.y_ - data.y_bar_) ;
  
  // Create the corresponding scaled/transformed data 
  RegressionData<mat> scaled_data(
      X_tilde,
      y_tilde,
      sp_mat(data.p_, data.p_),
      ones(data.p_), false, false) ;
  
  GroupLava<mat,MixedNorm::L1LINF> group_lava(scaled_data, Proj, group, regParam, control);
  
  List results = group_lava.solution_path(control);
  
  // Un-normalized vector of dense coefficients
  mat b = (C_inv * V * diagmat(D/(square(D) + 1)) * U.t()) * ( 
    (data.y_ - data.y_bar_) * ones(1, group_lava.lambdas_.size()) - 
      (data.X_.each_row() - data.X_bar_.t()) * group_lava.coefficients()
  ) ;

  group_lava.post_treatment(data, b) ;
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = group_lava.lambdas_,
      Named("l2") = gamma
    ),
    
    Named("coef")          = group_lava.coef_,
    Named("coef_debiased") = diagmat(1/data.norm_X_) * b + group_lava.debiased_coefficients(),
    Named("sparse_coef")   = group_lava.sparse_coef_,
    Named("sparse_coef_debiased") = group_lava.debiased_coefficients(),
    Named("active")        = group_lava.active_var(),
    Named("intercept")     = group_lava.intercept_,
    Named("intercept_debiased") = group_lava.intercept_debiased_,
    Named("normx")         = group_lava.data_.norm_X_,
    Named("df")            = group_lava.df_,
    Named("monitoring")    = results
  
  );
  
}

// [[Rcpp::export]]
List group_lava_coop_dense_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // Boolean for intercept
    const arma::uvec  &group    , // Vector of groups
    const List        &regParam    , // regularization parameters
    const List        &control    // config of the optimization 
) {
  
  // Data declaration and standardization
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
  vec lambda_factor = as<vec>(regParam["lambda_factor"]) ;
  double gamma = as<double>(regParam["gamma"]) ;
  data.scale_struct(sqrt(gamma)*ones(data.p_)) ;
  // data.scale_struct(sqrt(gamma)*pow(lambda_factor,-1)) ;
  // data.scale_regressors(lambda_factor) ;
  
  // Create the scaled/transformed data
  mat C_inv = solve(trimatu(chol(data.S_.as_dense())), eye(data.p_, data.p_)) ;
  mat U, V ; vec D ; // U D V = X C^-1
  svd_econ(U, D, V, (data.X_.each_row() - data.X_bar_.t())*C_inv) ;
  mat Proj = U * diagmat(square(D)/(square(D) + 1)) * U.t() ;
  mat K12 = U * diagmat(1/sqrt(square(D) + 1)) * U.t() ;
  mat X_tilde = K12 * (data.X_.each_row() - data.X_bar_.t()) ;
  vec y_tilde = K12 * (data.y_ - data.y_bar_) ;
  
  // Create the corresponding scaled/transformed data 
  RegressionData<mat> scaled_data(
      X_tilde,
      y_tilde,
      sp_mat(data.p_, data.p_),
      ones(data.p_), false, false) ;
  
  GroupLava<mat,MixedNorm::COOP> group_lava(scaled_data, Proj, group, regParam, control);
  
  List results = group_lava.solution_path(control);
  
  // Un-normalized vector of dense coefficients
  mat b = (C_inv * V * diagmat(D/(square(D) + 1)) * U.t()) * ( 
    (data.y_ - data.y_bar_) * ones(1, group_lava.lambdas_.size()) - 
      (data.X_.each_row() - data.X_bar_.t()) * group_lava.coefficients()
  ) ;
  
  group_lava.post_treatment(data, b) ;
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = group_lava.lambdas_,
      Named("l2") = gamma
    ),
    
    Named("coef")          = group_lava.coef_,
    Named("coef_debiased") = diagmat(1/data.norm_X_) * b + group_lava.debiased_coefficients(),
    Named("sparse_coef")   = group_lava.sparse_coef_,
    Named("sparse_coef_debiased") = group_lava.debiased_coefficients(),
    Named("active")        = group_lava.active_var(),
    Named("intercept")     = group_lava.intercept_,
    Named("intercept_debiased") = group_lava.intercept_debiased_,
    Named("normx")         = group_lava.data_.norm_X_,
    Named("df")            = group_lava.df_,
    Named("monitoring")    = results
  
  );
  
}
