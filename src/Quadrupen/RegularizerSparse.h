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
#include "PenaltySimple.h"
#include "PenaltyGroup.h"

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
  
  vec   nzeros_ ; // contains non-zero values of all betas (for all lambda values)
  vec debiased_ ; // contains debiased non-zero values of all betas (for all lambda values)
  vector<uvec> active_ ; // successively activated variable (for all lambda values)
  vector<double >intercept_debiased_ ; // debiased vector of intercept values (for all lambda values)
  vec beta_debiased_ ; // vector of current active beta debiased (for the current lambda)
  
  const sp_mat coefficients() const {
    vector<uword> rowA, colA ;
    uword current_col = 0;
    for (const auto& a : active_) {
      rowA.insert(rowA.end(), a.begin(), a.end()) ;
      colA.insert(colA.end(), a.n_elem, current_col) ;
      current_col++;
    }
    return sp_mat(join_cols(urowvec(rowA), urowvec(colA)),
                  nzeros_, data_.p_, active_.size(), true, false) ;
  }
  
  const sp_mat debiased_coefficients() const { 
    vector<uword> rowA, colA ;
    uword current_col = 0;
    for (const auto& a : active_) {
      rowA.insert(rowA.end(), a.begin(), a.end()) ;
      colA.insert(colA.end(), a.n_elem, current_col) ;
      current_col++;
    }
    return sp_mat(join_cols(urowvec(rowA), urowvec(colA)),
                  debiased_, data_.p_, active_.size(), true, false) ;
  }
  
  const sp_mat active_var() const { 
    vector<uword> rowA, colA ;
    uword current_col = 0;
    for (const auto& a : active_) {
      rowA.insert(rowA.end(), a.begin(), a.end()) ;
      colA.insert(colA.end(), a.n_elem, current_col) ;
      current_col++;
    }
    return sp_mat(join_cols(urowvec(rowA), urowvec(colA)),
                  vec(rowA.size(), fill::ones), data_.p_, active_.size(), true, false) ;
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
  using Regularizer<matrix>::lambda_factor_ ;
  
  SimpleSparseRegularizer() {} ;
  SimpleSparseRegularizer(const RegressionData<matrix>&, const List&);
  double get_lambda_max() 
  {
    return(penalty_.dual_norm(data_.XTy_, lambda_factor_));
  }
  
  ActiveSet<matrix> set_       ; // Active set of variable and data
  SimplePenalty<norm> penalty_ ; // main penalty object 
  
};

template <typename matrix, SimpleNorm norm>
SimpleSparseRegularizer<matrix, norm>::SimpleSparseRegularizer(
    const RegressionData<matrix>& data, const List& regParam) : 
  SparseRegularizer<matrix>::SparseRegularizer(data, regParam) {}

// ====================================================
// Group-Sparse Regularizers       

template <typename matrix, GroupNorm norm>
class GroupSparseRegularizer : public SparseRegularizer<matrix> {
public:
  
  using Regularizer<matrix>::data_ ;
  using Regularizer<matrix>::lambda_factor_ ;
  
  GroupSparseRegularizer() {} ;
  GroupSparseRegularizer(const RegressionData<matrix>&, const List&);
  double get_lambda_max() 
  {
    return(
      penalty_.dual_norm(data_.XTy_, set_.grp_sizes_, lambda_factor_)
    );
  }
  
  ActiveSetGroup<matrix> set_ ; // Active set of variable and data
  GroupPenalty<norm> penalty_ ; // main penalty object 

};

template <typename matrix, GroupNorm norm>
GroupSparseRegularizer<matrix, norm>::GroupSparseRegularizer(
    const RegressionData<matrix>& data, const List& regParam) : 
  SparseRegularizer<matrix>::SparseRegularizer(data, regParam) {}

#endif

