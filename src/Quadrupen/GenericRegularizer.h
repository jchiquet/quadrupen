/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _GenericRegularizer_H
#define _GenericRegularizer_H

#include "RegressionData.h"
#include "Penalty.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, Norm norm>
class GenericRegularizer {
public:
  
  GenericRegularizer() {} ;
  GenericRegularizer(const RegressionData<matrix>&, const List&);
  
  double get_lambda_max() {
    return(penalty_.dual_norm(data_.XTy_));
  } ;
  
  void get_lambda_seq(const List&);
  
  // Getter functions to access private members
  const RegressionData<matrix>& data()    const { return data_      ; }
  const matrix& coefficients()            const { return coef_      ; }
  const vector<double>& intercept()       const { return intercept_ ; }
  const vector<double>& path_tuning()     const { return lambdas_   ; }
  const vector<double>& degrees_freedom() const { return df_        ; }
  
  RegressionData<matrix> data_ ; // data structure
  Penalty<norm> penalty_       ; // main penalty object 
  vector<double> lambdas_      ; // vector of parameters tuning the man penalty
  vec lambda_factor_           ; // main penalty factors
  vector<double> intercept_    ; // vector of intercept term
  matrix coef_                 ; // matrix of coefficients
  vector<double> df_           ; // degrees of freedom along the path
  uvec all                     ; // a vector with all variable indices
};

template <typename matrix, Norm norm>
GenericRegularizer<matrix, norm>::GenericRegularizer(
  const RegressionData<matrix>& data, const List& regParam) : data_ (data)
{
  all = regspace<uvec>(0,data_.p_-1) ;
}

template <typename matrix, Norm norm>
void GenericRegularizer<matrix,norm>::get_lambda_seq(const List& regParam) {
  if (regParam[0] != R_NilValue) {
    lambdas_  = as<vector<double>>(regParam["lambda"]) ;
  } else {
    double lambda_max = this->get_lambda_max() ;
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

