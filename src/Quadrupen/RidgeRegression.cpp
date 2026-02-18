/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "RidgeRegression.h"

using namespace Rcpp;
using namespace arma;

RidgeRegression::RidgeRegression(
  RegressionData<mat>& data, const bool& intercept, const List& regParam) :
  GenericRegularizer::GenericRegularizer(data, intercept, regParam) 
{
  penalty_ = Penalty<Norm::RIDGE>() ;
  lambda_factor_ = as<vec>(regParam["lambda_factor"]) ;
  get_lambda_seq(regParam) ;
}

List RidgeRegression::solution_path(const mat& C_inv) {

  // SVD DECOMPOSITION OF ( X * C^-1)
  vec eta ; // eigen value of X cinv
  mat U   ; // left singular vectors of X
  mat V   ; // right singular vectors of X
  svd_econ(U, eta, V, (data_.X_.each_row() - data_.X_bar_.t())*C_inv) ;
  
  mat C_invV = C_inv * V ;
  mat Uty = trans(U) * data_.y_ ;

  vector<double> timing ; // successive timing for solving for each lambda value
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda : lambdas_) {
    // computing the structured ridge estimate
    vec beta = (C_invV * diagmat(eta/(square(eta) + lambda)) * Uty) / data_.norm_X_ ;
    coef_.push_back(beta);
    // estimating the intercept term
    const_.push_back(data_.y_bar_ - dot(beta, data_.X_bar_));  
    // computing the estimated degrees of freedom
    df_.push_back(sum(square(eta)/(square(eta) + lambda)));
    
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