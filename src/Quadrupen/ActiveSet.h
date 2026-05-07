/*
 * Author: Julien CHIQUET
 * MIA PS
 */
#pragma once

using namespace Rcpp;
using namespace arma;

#include "RegressionData.h"

template <typename matrix> 
class ActiveSet {

public:
  
  // VARIABLES FOR HANDLING THE ACTIVE SET
  uvec A_           ; // set of currently activated variables
  uvec is_in_       ; // indicator of active variables (0/1)
  mat XATXA_, XTXA_ ; // matrices of currently activated variables
  mat XATXAinv_     ; 
  bool use_chol_    ; // Maintain a Cholesky factorization along the active set algorithm
  mat R_            ; // Cholesky decomposition of XATXA
  
  ActiveSet() {} ;
  ActiveSet(const RegressionData<matrix> &data, const bool use_chol=true) ;
  ActiveSet(const RegressionData<matrix> &data, const uvec&, const bool use_chol) ;
  
  // ACTIVE SET HANDLING
  void add_var(uword, const RegressionData<matrix> &) ; // add a variable in the active set
  void add_vars(uvec, const RegressionData<matrix> &) ; // add a list of variables in the active set
  void del_var(uword, vec&) ; // remove the variable activated in position ind_var_out
  void del_vars(uvec, vec&) ; // remove a set of non contiguous variables
  void reset() ; // empty the active set 
  const uword size() const { return A_.n_elem ; }
  
  // Update Cholesky factorisation by inserting the last activated variables
  void update_Cholesky() ; 
  // Update Cholesky factorisation by inserting the last n_new activated variables
  void update_Cholesky_block(uword n_new) ;
                                                  
  // Downdate Cholesky factorisation by removing the specified variables
  void downdate_Cholesky(uword j) ; 
  
  // Inverse the currently active Gram matrix
  void inverse_Gram() ; 
  
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
  uword n_new = vars.n_elem;
  uword p_old = size();
  
  for(uword v : vars) is_in_[v] = 1;
  A_.resize(p_old + n_new);
  A_.tail(n_new) = vars;

  // new columns for XTXA_ (p x n_new)
  mat new_cols = data.X_.t() * data.X_.cols(vars) - 
    data.n_ * data.X_bar_ * data.X_bar_.rows(vars).t() + 
    data.S_.cols(vars);

  if (p_old == 0) { // If empty matrix
    XTXA_  = new_cols;
    XATXA_ = new_cols.rows(vars);
  } else { // Add to existing set of variables
    mat bottom_left = XTXA_.rows(vars); 
    
    XTXA_ = join_rows(XTXA_, new_cols);
    // [ XATXA_old   | bottom_left.t() ]
    // [ bottom_left | new_cols.rows(vars) ]
    mat bottom_right = new_cols.rows(vars);
    XATXA_ = join_cols(join_rows(XATXA_, bottom_left.t()), 
                       join_rows(bottom_left, bottom_right));
  }
  if (use_chol_) update_Cholesky_block(n_new) ;
}

template <typename matrix>
void ActiveSet<matrix>::del_var(uword ivar_out, vec& beta) {
  is_in_[A_[ivar_out]] = 0  ; // update the active set
  A_.shed_row(ivar_out)     ;
  XTXA_.shed_col(ivar_out)  ;
  XATXA_.shed_col(ivar_out) ;
  XATXA_.shed_row(ivar_out) ;
  beta.shed_row(ivar_out)   ;
  
  if (use_chol_) downdate_Cholesky(ivar_out) ;
  
}

template <typename matrix>
void ActiveSet<matrix>::del_vars(uvec ivars, vec& beta) {
  ivars = sort(ivars, "descend");
  for (uword i=0 ; i <ivars.n_elem ; i++) {
    del_var(ivars[i], beta) ;
  }
}

template <typename matrix>
void ActiveSet<matrix>::update_Cholesky() {
  uword p = XATXA_.n_cols ;
  
  if (p == 1) {
    R_ = sqrt(XATXA_);
  } else {
    colvec rp  = zeros<colvec>(p,1);
    rp.subvec(0, p-2) = solve(trimatu(R_.submat(0,0,p-2,p-2)).t(), 
              XATXA_.submat(0, p-1, p-2, p-1), 
              solve_opts::fast);
    rp(p-1) = sqrt(XATXA_(p-1,p-1) - dot(rp,rp));
    R_ = join_rows( join_cols(R_, zeros<mat>(1,p-1)) , rp);
  }
}

template <typename matrix>
void ActiveSet<matrix>::update_Cholesky_block(uword n_new) {
  uword p_total = XATXA_.n_cols;
  uword p_old = p_total - n_new;
  
  if (p_old == 0) {
    R_ = chol(XATXA_);
  } else {
    // Solve R_old.t() * R_new_block = XATXA_joint
    // XATXA_joint is the upper right block (p_old x n_new)
    mat R_new_cols = solve(trimatu(R_).t(), 
                           XATXA_.submat(0, p_old, p_old - 1, p_total - 1), 
                           solve_opts::fast);
    
    // Compute the Schur complement for the new diagonal block
    // S = XATXA_new - R_new_cols.t() * R_new_cols
    mat S = XATXA_.submat(p_old, p_old, p_total - 1, p_total - 1) - 
      R_new_cols.t() * R_new_cols;
    mat R_bottom_right = chol(S);
    
    // Final assembly
    // [ R_old | R_new_cols     ]
    // [ 0     | R_bottom_right ]
    R_ = join_rows(join_cols(R_, zeros<mat>(n_new, p_old)), 
                   join_cols(R_new_cols, R_bottom_right));
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
      r = std::hypot(x(0), x(1));
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

template <typename matrix>
void ActiveSet<matrix>::inverse_Gram() {
  if (use_chol_) {
    XATXAinv_ = solve(trimatu(R_), solve(trimatl(R_.t()), eye<mat>(R_.n_cols, R_.n_cols)));
  } else {
    XATXAinv_ = inv_sympd(XATXA_, inv_opts::allow_approx);
  }
}



