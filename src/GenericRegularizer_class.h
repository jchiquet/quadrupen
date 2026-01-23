/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _GenericRegularizer_H
#define _GenericRegularizer_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "RegressionData_class.h"
#include "Penalty_class.h"

#define ZERO 2e-16

using namespace Rcpp;
using namespace arma;
using namespace std;

class GenericRegularizer {
public:
  GenericRegularizer(RegressionData&, const bool&, const List&);
  virtual ~GenericRegularizer() = 0 ;
  
  double get_lambda_max() {
    return(penalty_.dual_norm(data_.XTy_));
  } ;
  
  void lambda_seq(const List&);
  
  // Getter functions to access private members
  const RegressionData& data()            const { return data_      ; }
  const bool& intercept()                 const { return intercept_ ; }
  const vector<vec>& coefficients()       const { return beta_      ; }
  const vector<double>& interceptTerm()   const { return mu_        ; }
  const vector<double>& tuning_param()    const { return lambda_    ; }
  const vector<double>& degrees_freedom() const { return df_        ; }
  
protected:
  
  // Variables common to all regularizers
  RegressionData data_    ; // data structure
  bool intercept_         ; // does the model include intercept
  Penalty penalty_        ; // main penalty object 
  vector<double> lambda_  ; // vector of parameters tuning the man penalty
  vec lambda_factor_      ; // main penalty factors
  vector<vec> beta_       ; // matrix of coefficients
  vector<double> df_      ; // degrees of freedom along the path
  vector<double> mu_      ; // vector of intercept term

};

#endif

