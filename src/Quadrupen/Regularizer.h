/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _Regularizer_H
#define _Regularizer_H

#include "RegressionData.h"
#include "Penalty.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class Regularizer {
public:
  
  Regularizer() {} ;
  Regularizer(const RegressionData<matrix>&, const List&);
  
  void get_lambda_seq(double, const List&);
  
  // Getter functions to access private members
  const RegressionData<matrix>& data()    const { return data_      ; }
  const matrix& coefficients()            const { return coef_      ; }
  const vector<double>& intercept()       const { return intercept_ ; }
  const vector<double>& path_tuning()     const { return lambdas_   ; }
  const vector<double>& degrees_freedom() const { return df_        ; }
  
  RegressionData<matrix> data_ ; // data structure
  vector<double> lambdas_      ; // vector of parameters tuning the man penalty
  vec lambda_factor_           ; // main penalty factors
  vector<double> intercept_    ; // vector of intercept term
  matrix coef_                 ; // matrix of coefficients
  vector<double> df_           ; // degrees of freedom along the path
  uvec all                     ; // a vector with all variable indices
};

template <typename matrix>
Regularizer<matrix>::Regularizer(
  const RegressionData<matrix>& data, const List& regParam) : data_ (data)
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

template <typename matrix, SimpleNorm norm>
class SimpleRegularizer : public Regularizer<matrix> {
public:
  
  SimpleRegularizer() {} ;
  SimpleRegularizer(const RegressionData<matrix>&, const List&);
  
  SimplePenalty<norm> penalty_ ; // main penalty object 
};

template <typename matrix, SimpleNorm norm>
SimpleRegularizer<matrix, norm>::SimpleRegularizer(
    const RegressionData<matrix>& data, const List& regParam) : 
  Regularizer<matrix>::Regularizer(data, regParam) {}


template <typename matrix, MixedNorm norm>
class GroupRegularizer : public Regularizer<matrix> {
public:
  
  GroupRegularizer() {} ;
  GroupRegularizer(const RegressionData<matrix>&, const List&);
  
  MixedPenalty<norm> penalty_ ; // main penalty object 
};

template <typename matrix, MixedNorm norm>
GroupRegularizer<matrix, norm>::GroupRegularizer(
    const RegressionData<matrix>& data, const List& regParam) : 
  Regularizer<matrix>::Regularizer(data, regParam) {}


#endif

