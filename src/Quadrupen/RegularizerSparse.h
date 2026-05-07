/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "Regularizer.h"
#include "RegressionData.h"

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


