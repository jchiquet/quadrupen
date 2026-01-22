/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_RIDGE_H
#define _quadrupen_RIDGE_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "RegressionData_class.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

class RidgeRegression {
 public:
   RidgeRegression(RegressionData&, const bool&, const List&);

  void lambda_seq(const List&);

  List solution_path(const mat&);

  // getter functions to access private members
  const RegressionData& data()            const { return data_      ; }
  const bool& intercept()                 const { return intercept_ ; }
  const vector<vec>& coefficients()       const { return beta_      ; }
  const vector<double>& interceptTerm()   const { return mu_        ; }
  const vector<double>& tuning_param()    const { return lambda_    ; }
  const vector<double>& degrees_freedom() const { return df_        ; }

 private:
   RegressionData data_    ; // data structure
   bool intercept_         ; // does the model include intercept
   vector<double> lambda_  ; // vector of ridge tuning parameters
   vector<vec> beta_       ; // matrix of coefficients
   vector<double> df_      ; // degrees of freedom along the path
   vector<double> mu_      ; // vector of intercept term

};

#endif

