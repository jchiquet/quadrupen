/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/GroupLava.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List group_lava_l1l2_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {

  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  GroupLava<mat,GroupSparseNorm::L1L2> group_lava(data, group, regParam, control) ;
  List results = group_lava.solution_path(control) ;
  group_lava.post_treatment() ;

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = group_lava.lambdas_,
      Named("l2") = as<double>(regParam["gamma"])
    ),
    Named("coef")                 = group_lava.coef_,
    Named("coef_debiased")        = diagmat(1/data.norm_X_) * group_lava.b_ + group_lava.debiased_coefficients(),
    Named("sparse_coef")          = group_lava.sparse_coef_,
    Named("sparse_coef_debiased") = group_lava.debiased_coefficients(),
    Named("active")               = group_lava.active_var(),
    Named("intercept")            = group_lava.intercept_,
    Named("intercept_debiased")   = group_lava.intercept_debiased_,
    Named("normx")                = group_lava.data_.norm_X_,
    Named("df")                   = group_lava.df_,
    Named("monitoring")           = results
  ) ;

}

// [[Rcpp::export]]
List group_lava_l1linf_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {

  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  GroupLava<mat,GroupSparseNorm::L1LINF> group_lava(data, group, regParam, control) ;
  List results = group_lava.solution_path(control) ;
  group_lava.post_treatment() ;

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = group_lava.lambdas_,
      Named("l2") = as<double>(regParam["gamma"])
    ),
    Named("coef")                 = group_lava.coef_,
    Named("coef_debiased")        = diagmat(1/data.norm_X_) * group_lava.b_ + group_lava.debiased_coefficients(),
    Named("sparse_coef")          = group_lava.sparse_coef_,
    Named("sparse_coef_debiased") = group_lava.debiased_coefficients(),
    Named("active")               = group_lava.active_var(),
    Named("intercept")            = group_lava.intercept_,
    Named("intercept_debiased")   = group_lava.intercept_debiased_,
    Named("normx")                = group_lava.data_.norm_X_,
    Named("df")                   = group_lava.df_,
    Named("monitoring")           = results
  ) ;

}

// [[Rcpp::export]]
List group_lava_coop_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {

  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  GroupLava<mat,GroupSparseNorm::COOP> group_lava(data, group, regParam, control) ;
  List results = group_lava.solution_path(control) ;
  group_lava.post_treatment() ;

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = group_lava.lambdas_,
      Named("l2") = as<double>(regParam["gamma"])
    ),
    Named("coef")                 = group_lava.coef_,
    Named("coef_debiased")        = diagmat(1/data.norm_X_) * group_lava.b_ + group_lava.debiased_coefficients(),
    Named("sparse_coef")          = group_lava.sparse_coef_,
    Named("sparse_coef_debiased") = group_lava.debiased_coefficients(),
    Named("active")               = group_lava.active_var(),
    Named("intercept")            = group_lava.intercept_,
    Named("intercept_debiased")   = group_lava.intercept_debiased_,
    Named("normx")                = group_lava.data_.norm_X_,
    Named("df")                   = group_lava.df_,
    Named("monitoring")           = results
  ) ;

}
