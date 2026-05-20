/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "RegressionData.h"

using arma::vec;
using arma::uword;
using arma::conv_to;
using arma::logspace;
using Rcpp::List;
using Rcpp::as;
using std::vector;

template <typename matrix>
class Regularizer {
public:
  
  Regularizer() {} ;
  Regularizer(const RegressionData<matrix>&, const List&);

  RegressionData<matrix> data_ ; // data structure
  double gamma_                ; // overall amount of minor penalty (not leading the path)
  vector<double> lambdas_      ; // vector of parameters tuning the main penalty
  vec lambda_factor_           ; // weights for the main penalty
  vector<double> intercept_    ; // vector of intercept term
  matrix coef_                 ; // matrix of coefficients
  vec beta_                    ; // vector of current parameters (for fix lambda value)
  vec grad_                    ; // vector of current gradient (smooth part)
  vector<double> df_           ; // degrees of freedom along the path

  void get_lambda_seq(double, const List&);
  
  // Getter functions to access private members
  const RegressionData<matrix>& data()    const { return data_      ; }
  const matrix& coefficients()            const { return coef_      ; }
  const vector<double>& intercept()       const { return intercept_ ; }
  const vector<double>& path_tuning()     const { return lambdas_   ; }
  const vector<double>& degrees_freedom() const { return df_        ; }
};

template <typename matrix>
Regularizer<matrix>::Regularizer(
  const RegressionData<matrix>& data, const List& regParam) : 
  data_ (data),
  gamma_(as<double>(regParam["gamma"])),
  lambda_factor_(as<vec>(regParam["lambda_factor"]))
{}

template <typename matrix>
void Regularizer<matrix>::get_lambda_seq(double lambda_max, const List& regParam) {
  if (regParam["lambda"] != R_NilValue) {
    lambdas_  = as<vector<double>>(regParam["lambda"]) ;
  } else {
    lambdas_ = conv_to<vector<double>>::from(
      logspace(
        log10(lambda_max),
        log10(as<double>(regParam["min_ratio"])*lambda_max),
        as<uword>(regParam["n_lambda"])
      )
    );
  }
}

