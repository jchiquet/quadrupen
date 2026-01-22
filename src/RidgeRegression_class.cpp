/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "RidgeRegression_class.h"

using namespace Rcpp;
using namespace arma;

RidgeRegression::RidgeRegression(
  RegressionData& data, const bool& intercept, const List& regParam) :
  data_ (data), intercept_ (intercept)
{
  lambda_seq(regParam) ;
}


void RidgeRegression::lambda_seq(const List& regParam) {
  if (regParam["l2"] != R_NilValue) {
    lambda_  = as<vector<double>>(regParam["l2"]) ;
  } else {
    double lambda_max = as<double>(regParam["lambda_max"]);
    lambda_ = conv_to<vector<double>>::from(
      logspace(log10(lambda_max),
      log10(as<double>(regParam["min_ratio"])*lambda_max),
      as<uword>(regParam["n_lambda"])
    ));
  }
}

List RidgeRegression::solution_path(const mat& C_inv) {

  // SVD DECOMPOSITION OF ( X * C^-1)
  arma::vec eta ; // eigen value of X cinv
  arma::mat U   ; // left singular vectors of X
  arma::mat V   ; // right singular vectors of X
  svd_econ(U, eta, V, (data_.X_.each_row() - data_.X_bar_.t())*C_inv) ;

  arma::mat C_invV = C_inv * V ;
  arma::mat Uty = trans(U) * data_.y_ ;

  vector<double> timing ; // successive timing for solving for each lambda value
  wall_clock timer ; timer.tic(); // clock
  for(auto current_lambda : lambda_) {
    // computing the structured ridge estimate
    vec current_beta = (C_invV * diagmat(eta/(square(eta) + current_lambda)) * Uty) / data_.norm_X() ;
    beta_.push_back(current_beta);
    // estimating the intercept term
    mu_.push_back(as_scalar(data_.y_bar_ - dot(current_beta, data_.X_bar_)));  
    // computing the estimated degrees of freedom
    df_.push_back(sum(square(eta)/(square(eta) + current_lambda)));
    
    timing.push_back(timer.toc()) ;
  }

  return(
    List::create(
      Named("max_grd")        = R_NaReal,
      Named("convergence")    = 0,
      Named("pensteps_timer") = timing 
    )
  );
  
}