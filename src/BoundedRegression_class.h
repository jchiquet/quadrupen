/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _BoundedRegression_H
#define _BoundedRegression_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "RegressionData_class.h"
#include "Penalty_class.h"
#include "Optimizer_class.h"
#include "utils.h"

#define ZERO 2e-16 // practical zero

using namespace Rcpp;
using namespace arma;
using namespace std;

class BoundedRegression {
public:
  BoundedRegression(RegressionData&, const bool&, const List&);
  
  void lambda_seq(const List&);
  
  List solution_path(const List&);
  
  // Getter functions to access private members
  const RegressionData& data()            const { return data_      ; }
  const bool& intercept()                 const { return intercept_ ; }
  const vector<vec>& coefficients()       const { return beta_      ; }
  const vector<double>& interceptTerm()   const { return mu_        ; }
  const vector<double>& tuning_param()    const { return lambda_    ; }
  const vector<double>& degrees_freedom() const { return df_        ; }
  
private:

  // Variables common to all regularizers
  RegressionData data_    ; // data structure
  Penalty penalty_        ; // Penalty object 
  vector<double> lambda_  ; // vector of penalty tuning parameters
  bool intercept_         ; // does the model include intercept
  vector<vec> beta_       ; // matrix of coefficients
  vector<double> df_      ; // degrees of freedom along the path
  vector<double> mu_      ; // vector of intercept term
  
  // Specific to Bounded regression
  uvec all          ; // remains constant
  mat    XTX        ; // Gram matrix
  vec linf_weights  ; // l_inf penalty factors
  double l2_penalty ; // overall amount of l2 penalty
  sp_mat S_scaled   ; // structuring matrix scaled according to l_inf weights
  uvec B            ; // guys reaching the boundary
  urowvec iB_       ; // contains row indices of the bounded variables
  urowvec jB_       ; // contains column indices of the bounded variables
  
  // Compute degrees of freedom for the current estimate
  double get_df() ; 
};

#endif

