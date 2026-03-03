/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _Regularizer_H
#define _Regularizer_H

#include "RegressionData.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class Regularizer {
public:
  
  Regularizer() {} ;
  Regularizer(const RegressionData<matrix>&, const List&);

  RegressionData<matrix> data_ ; // data structure
  double gamma_                ; // overall amount of minor penalty (not leading the path)
  vector<double> lambdas_      ; // vector of parameters tuning the main penalty
  vec lambda_factor_           ; // main penalty factors
  vector<double> intercept_    ; // vector of intercept term
  matrix coef_                 ; // matrix of coefficients
  vec beta_                    ; // vector of current parameters
  vec grad_                    ; // vector of current gradient (smooth part)
  vector<double> df_           ; // degrees of freedom along the path
  uvec all                     ; // a vector with all variable indices
  
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
  data_ (data), gamma_(as<double>(regParam["gamma"])), lambda_factor_(as<vec>(regParam["lambda_factor"]))
{
  all = regspace<uvec>(0,data_.p_-1) ;
}

template <typename matrix>
void Regularizer<matrix>::get_lambda_seq(double lambda_max, const List& regParam) {
  if (regParam[0] != R_NilValue) {
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

#endif

