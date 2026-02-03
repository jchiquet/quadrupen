/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _ElasticNet_H
#define _ElasticNet_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "GenericRegularizer_class.h"
#include "ActiveSet_class.h"
#include "Optimizer_class.h"

#define ZERO 2e-16 // practical zero

using namespace Rcpp;
using namespace arma;
using namespace std;

class ElasticNet : public GenericRegularizer<mat> {

  public:
  ElasticNet(RegressionData<mat>&, const bool&, const List&, const List&);
  
  List solution_path(const List&);
  
  const double& struct_tuning() const { return gamma_ ; }

  const sp_mat coefficients() const { 
    return sp_mat(join_cols(iA_, jA_), nzeros_, lambdas_.size(), data_.p_) ; 
  }

  const sp_mat active_var() const { 
    return sp_mat(join_cols(iA_, jA_), vec(iA_.n_elem, fill::ones), lambdas_.size(), data_.p_) ; 
  }
  
private:

  // Specific to Elastic-Net regularization
  double gamma_ ; // overall amount of l2 penalty
  vec beta_     ; // vector of current parameters
  vec grad_     ; // vector of current gradient (smooth part)
  vec   nzeros_ ; // contains non-zero value of beta
  urowvec iA_   ; // contains row indices of the non-zero values
  urowvec jA_   ; // contains column indices of the non-zero values

  // Compute degrees of freedom for the current estimate
  double get_df() ; 

};

#endif

