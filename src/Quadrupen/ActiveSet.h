/*
 * Author: Julien CHIQUET
 */
#ifndef _ActiveSet_H
#define _ActiveSet_H

using namespace Rcpp;
using namespace arma;

#include "RegressionData.h"

template <typename matrix> class ActiveSet {

public:
  
  // VARIABLES FOR HANDLING THE ACTIVE SET
  uvec A_           ; // set of currently activated variables
  uvec is_in_       ; // indicator of active variables (0/1)
  mat XATXA_, XTXA_ ; // matrices of currently activated variables
  bool use_chol_    ; // Maintain a Cholesky factorization along the active set algorithm
  mat R_            ; // Cholesky decomposition of XATXA
  
  ActiveSet() {} ;
  ActiveSet(const RegressionData<matrix> &data, const bool use_chol=true) ;
  ActiveSet(const RegressionData<matrix> &data, const uvec&, const bool use_chol) ;
  
  // ACTIVE SET HANDLING
  void add_var(uword, const RegressionData<matrix> &) ; // add a variable in the active set
  void add_vars(uvec, const RegressionData<matrix> &) ; // add a list of variables in the active set
  void del_var(uword) ; // remove the variable activated in position ind_var_out
  void del_vars(uvec) ; // remove a set of non contiguous varables
  void reset()        ; // empty the active set 
  const uword size() const { return A_.n_elem ; }
  
  // Update Cholesky factorisation by inserting the last activated variables
  void update_Cholesky() ; 
  
  // Downdate Cholesky factorisation by removing the specified variables
  void downdate_Cholesky(uword j) ; 
};

template <typename matrix>
ActiveSet<matrix>::ActiveSet(const RegressionData<matrix>& data, const bool use_chol) :
  use_chol_(use_chol) {
  is_in_.zeros(data.p_) ;
}

template <typename matrix>
ActiveSet<matrix>::ActiveSet(const RegressionData<matrix>& data, const uvec& A0, const bool use_chol) :
  use_chol_(use_chol) {
  is_in_.zeros(data.p_) ;
  add_vars(A0, data)    ;
}

template <typename matrix>
void ActiveSet<matrix>::reset() {
  A_.reset()      ;
  is_in_.zeros()  ;
  XATXA_.reset()  ;
  XTXA_.reset()   ; 
  R_.reset()      ;
}

template <typename matrix>
void ActiveSet<matrix>::add_var(uword var_in, const RegressionData<matrix>& data) {
  A_.resize(size() + 1) ;
  A_.tail(1) = var_in   ;
  is_in_[var_in] = 1    ;

  vec new_col = data.X_.t() * data.X_.col(var_in) - 
    data.n_ * data.X_bar_ * as_scalar(data.X_bar_[var_in]) + data.S_.col(var_in) ;

  if (!XATXA_.is_empty()) {
    XATXA_ = join_cols(XATXA_, XTXA_.row(var_in)) ;
  }
  XTXA_  = join_rows(XTXA_ , new_col) ;
  XATXA_ = join_rows(XATXA_, trans(XTXA_.row(var_in))) ;

  if (use_chol_) update_Cholesky() ;
}

template <typename matrix>
void ActiveSet<matrix>::add_vars(uvec vars, const RegressionData<matrix>& data) {
  for (uword v=0; v < vars.n_elem ; v++) {
    add_var(vars[v], data) ;
  }
}

template <typename matrix>
void ActiveSet<matrix>::del_var(uword ivar_out) {
  is_in_[A_[ivar_out]] = 0   ; // update the active set
  A_.shed_row(ivar_out)      ;
  XTXA_.shed_col(ivar_out)   ;
  XATXA_.shed_col(ivar_out)  ;
  XATXA_.shed_row(ivar_out)  ;

  if (use_chol_) downdate_Cholesky(ivar_out) ;
  
}

template <typename matrix>
void ActiveSet<matrix>::del_vars(uvec ivars) {
  ivars = sort(ivars, "descend");
  for (uword i=0 ; i <ivars.n_elem ; i++) {
    del_var(ivars[i]) ;
  }
}

template <typename matrix>
void ActiveSet<matrix>::update_Cholesky() {
  uword p = XATXA_.n_cols ;
  
  if (p == 1) {
    R_ = sqrt(XATXA_);
  } else {
    colvec rp  = zeros<colvec>(p,1);
    rp.subvec(0,p-2) = solve (trimatl(strans(R_)), XATXA_.submat(0,p-1,p-2,p-1));
    rp(p-1) = sqrt(XATXA_(p-1,p-1) - dot(rp,rp));
    R_ = join_rows( join_cols(R_, zeros<mat>(1,p-1)) , rp);
  }
  
}

template <typename matrix>
void ActiveSet<matrix>::downdate_Cholesky(uword j) {
  
  vec x = zeros<vec>(2,1);
  mat G = zeros<mat>(2,2);
  mat H = zeros<mat>(2,2);
  
  R_.shed_col(j);
  int p = R_.n_cols;
  double r;
  for (int k=j; k<p; k++) {
    x = R_.submat(k,k,k+1,k);
    
    if (x[1] != 0) {
      r = norm(x,2);
      G = {{x(0), x(1)}, {-x(1), x(0)}};
      G = G / r;
      x(0) = r; x(1) = 0;
    } else {
      G = eye(2,2);
    }
    R_.submat(k,k,k+1,k) = x;
    if (k < p-1) {
      R_.submat(k,k+1,k+1,p-1) = G * R_.submat(k,k+1,k+1,p-1);
    }
  }
  R_.shed_row(p);
}

#endif
