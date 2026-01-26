/*
 * Author: Julien CHIQUET
 */
#ifndef _ActiveSet_H
#define _ActiveSet_H

using namespace Rcpp;
using namespace arma;

#include "RegressionData_class.h"

template <typename matrix> class GenericRegularizer ;

template <typename matrix> class ActiveSet {
  friend class GenericRegularizer<matrix> ;
  friend class ElasticNet                 ;
  
protected:
  // VARIABLES FOR HANDLING THE ACTIVE SET
  vec  betaA_   ; // vector of currently activated features
  uvec A_       ; // set of currently activated variables
  uvec is_in_A_ ; // a vector to check if a variable is already in the active set
  mat XATXA_    ;
  mat XTXA_     ;
  vec grad_     ;
  
public:
  ActiveSet() {} ;
  ActiveSet(SEXP BETA0, RegressionData<matrix> &data) ;
  
  // ACTIVE SET HANDLING
  void add_var(uword, RegressionData<matrix> &) ; // add a variable in the active set
  void del_var(uword)  ; // remove the variable activated in position ind_var_out
  void add_vars(uword first_in, uword last_in, RegressionData<matrix> &)   ; // the same for a set 
  void del_vars(uword first_out, uword last_out) ; // of contiguous variables
  
  const vec& beta_active() const { return betaA_ ; }
  const uvec& active() const { return A_  ; }
  const uword& size() const { return A_.n_elem  ; }
  bool is_active (uword i) { return (is_in_A_(i)) ; }
  const mat& Gram_active() const { return XATXA_ ; }
  
};

template <typename matrix>
ActiveSet<matrix>::ActiveSet(SEXP BETA0, RegressionData<matrix>& data) {
  
  vec beta0 = as<vec>(BETA0)  ;
  A_    = find(beta0 != 0)    ;
  betaA_ = beta0.elem(A_)     ;
  is_in_A_.zeros(beta0.n_elem)  ;
  for (uword i=0; i<is_in_A_.n_elem;i++) {
    is_in_A_(A_(i)) = 1 ;
  }
  std::cout << "HERE" << std::endl;
  
  if (!betaA_.is_empty()) {
    XTXA_  = mat(data.X_.t() * data.X_.cols(A_)) - data.n_ * data.X_bar_ * (data.X_bar_.elem(A_)).t();
    // grad_  = -data_.XTy_ + XTXA_ * betaA_          ; // not sure grad should be here...
    XATXA_ = XTXA_.rows(A_)                        ;
  }
}

template <typename matrix>
void ActiveSet<matrix>::add_var(uword varIn, RegressionData<matrix>& data) {
  betaA_.resize(A_.n_elem+1); // update the vector of active parameters
  betaA_(A_.n_elem) = 0.0   ;
  A_.resize(A_.n_elem)      ; // update the active set
  A_[A_.n_elem] = varIn     ;
  is_in_A_[varIn] = 1       ;
  
  vec new_col = data.X_.t() * data.X_.col(varIn) - 
    data.n * data.X_bar_ * as_scalar(data.X_bar_[varIn]) ;
  // if (lambda2 > 0) {
  // Adding the column corresponding to the structurating matrix
  // new_col += spS.col(var_in);
  // }
  
  // UPDATE THE xtxA AND xAtxA MATRICES
  if (!XATXA_.is_empty()) {
    XATXA_ = join_cols(XATXA_, XTXA_.row(varIn)) ;
  }
  XTXA_  = join_rows(XTXA_ , new_col) ;
  XATXA_ = join_rows(XATXA_, trans(XTXA_.row(varIn))) ;

}

// must be sorted in increasing order
template <typename matrix>
void ActiveSet<matrix>::add_vars(uword first_varIn, uword last_varIn, RegressionData<matrix>& data) {
  for (uword k=first_varIn; k<=last_varIn; k++) {
    betaA_.resize(A_.n_elem+1); // update the vector of active parameters
    betaA_(A_.n_elem) = 0.0   ;
    A_.resize(A_.n_elem)      ; // update the active set
    A_[A_.n_elem] = k     ;
    is_in_A_[k] = 1       ;
  }
  
  mat new_cols = data.X.t() * data.X.cols(first_varIn, last_varIn) -
    data.n * data.X_bar_ * data.X_bar_.elem(first_varIn, last_varIn) ;

  // UPDATE THE xtxA AND xAtxA MATRICES
  if (!XATXA_.is_empty()) {
    XATXA_ = join_cols(XATXA_, XTXA_.rows(first_varIn, last_varIn)) ;
  }
  XTXA_  = join_rows(XTXA_, new_cols) ;
  XATXA_ = join_rows(XATXA_, trans(XTXA_.rows(first_varIn, last_varIn))) ;
}

template <typename matrix>
void ActiveSet<matrix>::del_var(uword indVarOut) {
  betaA_.shed_row(indVarOut)  ; // update the vector of active parameters
  is_in_A_[A_[indVarOut]] = 0 ; // update the active set
  A_.shed_row(indVarOut)      ;
  XTXA_.shed_col(indVarOut)   ;
  XATXA_.shed_col(indVarOut)  ;
  XATXA_.shed_row(indVarOut)  ;
}


template <typename matrix>
void ActiveSet<matrix>::del_vars(uword first_indVarOut, uword last_indVarOut) {
  for (uword k=last_indVarOut; k>=first_indVarOut; k--) {
    betaA_.shed_row(k)  ; // update the vector of active parameters
    is_in_A_[A_[k]] = 0 ; // update the active set
    A_.shed_row(k)      ;
    if (k == 0) {break;}
  }
  XTXA_.shed_cols(first_indVarOut, last_indVarOut)  ;
  XATXA_.shed_cols(first_indVarOut, last_indVarOut) ;
  XATXA_.shed_rows(first_indVarOut, last_indVarOut) ;
}

#endif
