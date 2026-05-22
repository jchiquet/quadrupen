/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/Lava.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List lava_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &control
) {

  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;

  Lava<mat,SparseNorm::L1> lava(data, regParam, control) ;
  List results = lava.solution_path(control) ;
  lava.post_treatment() ;

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = lava.lambdas_,
      Named("l2") = as<double>(regParam["gamma"])
    ),
    Named("coef")                 = lava.coef_,
    Named("coef_debiased")        = diagmat(1/data.norm_X_) * lava.b_ + lava.debiased_coefficients(),
    Named("sparse_coef")          = lava.sparse_coef_,
    Named("sparse_coef_debiased") = lava.debiased_coefficients(),
    Named("active")               = lava.active_var(),
    Named("intercept")            = lava.intercept_,
    Named("intercept_debiased")   = lava.intercept_debiased_,
    Named("normx")                = lava.data_.norm_X_,
    Named("df")                   = lava.df_,
    Named("monitoring")           = results
  ) ;

}
