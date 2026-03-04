/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 *         
 *  3 Classes for Sparse regularization:
 *  - SparseRegularizer
 *  - SimpleSparseRegularizer
 *  - GroupSparseRegularizer
 *  
 *  Inheritance:
 *  
 *  Regularizer-> SparseRegularizer -> SimpleSparseRegularizers
 *  Regularizer-> SparseRegularizer -> GroupSparseRegularizers
 */

#ifndef _RegularizerSparse_H
#define _RegularizerSparse_H

#include "Regularizer.h"
#include "RegressionData.h"
#include "ActiveSet.h"
#include "ActiveSetGroup.h"
#include "Penalty.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// ====================================================
// Sparse Regularizers 
// base class, common to SimpleSparse and GroupSparse

template <typename matrix>
class SparseRegularizer : public Regularizer<matrix> {
public:
  
  using Regularizer<matrix>::data_ ;
  using Regularizer<matrix>::lambdas_ ;
  
  SparseRegularizer() {} ;
  SparseRegularizer(const RegressionData<matrix>&, const List&);
  
  vec   nzeros_ ; // contains non-zero value of beta
  vec debiased_ ; // contains the debiased non-zero value of beta
  vector<double >intercept_debiased_ ; // contains the debiased vector of intercept
  urowvec iA_   ; // contains row indices of the non-zero values
  urowvec jA_   ; // contains column indices of the non-zero values
  
  const sp_mat coefficients() const {
    return sp_mat(join_cols(iA_, jA_), nzeros_, data_.p_, lambdas_.size()) ; 
  }
  
  const sp_mat debiased_coefficients() const { 
    return sp_mat(join_cols(iA_, jA_), debiased_, data_.p_, lambdas_.size()) ; 
  }
  
  const sp_mat active_var() const { 
    return sp_mat(join_cols(iA_, jA_), vec(iA_.n_elem, fill::ones), data_.p_, lambdas_.size()) ; 
  }

  const vector<double>& intercept_debiased() const { return intercept_debiased_ ; }

};

template <typename matrix>
SparseRegularizer<matrix>::SparseRegularizer(
  const RegressionData<matrix>& data, const List& regParam) : 
  Regularizer<matrix>(data, regParam)
  {}

// ====================================================
// Simple Sparse Regularizers

template <typename matrix, SimpleNorm norm>
class SimpleSparseRegularizer : public SparseRegularizer<matrix> {
public:
  
  using Regularizer<matrix>::data_ ;
  
  SimpleSparseRegularizer() {} ;
  SimpleSparseRegularizer(const RegressionData<matrix>&, const List&);
  double get_lambda_max() {return(penalty_.dual_norm(data_.XTy_));}
  
  ActiveSet<matrix> set_       ; // Active set of variable and data
  SimplePenalty<norm> penalty_ ; // main penalty object 
  
};

template <typename matrix, SimpleNorm norm>
SimpleSparseRegularizer<matrix, norm>::SimpleSparseRegularizer(
    const RegressionData<matrix>& data, const List& regParam) : 
  SparseRegularizer<matrix>::SparseRegularizer(data, regParam) {}

// ====================================================
// Group-Sparse Regularizers       

template <typename matrix, MixedNorm norm>
class GroupSparseRegularizer : public SparseRegularizer<matrix> {
public:
  
  using Regularizer<matrix>::data_ ;
  
  GroupSparseRegularizer() {} ;
  GroupSparseRegularizer(const RegressionData<matrix>&, const List&);
  double get_lambda_max() {return(penalty_.dual_norm(data_.XTy_, set_.grp_sizes_));}
  
  ActiveSetGroup<matrix> set_ ; // Active set of variable and data
  MixedPenalty<norm> penalty_ ; // main penalty object 

};

template <typename matrix, MixedNorm norm>
GroupSparseRegularizer<matrix, norm>::GroupSparseRegularizer(
    const RegressionData<matrix>& data, const List& regParam) : 
  SparseRegularizer<matrix>::SparseRegularizer(data, regParam) {}

#endif

