#include "active_set.hpp"

using namespace Rcpp;
using namespace arma;

// ______________________________________________________
// ACTIVE SET CLASS
//
// handle vector of active variable, with method for adding/removing
// elements

ACTIVE_SET::ACTIVE_SET(SEXP BETA0) {
  vec beta0 = as<vec>(BETA0)         ;
  activeVar = find(beta0 != 0)       ;
  betaActive = beta0.elem(activeVar) ;
  nbrIn = activeVar.n_elem           ;
  isActiveVar.zeros(beta0.n_elem)    ;
  for (int i=0; i<nbrIn;i++) {
    isActiveVar(activeVar(i)) = 1 ;
  }
}

void ACTIVE_SET::add_var(uword varIn) {
  betaActive.resize(nbrIn+1); // update the vector of active parameters
  betaActive(nbrIn) = 0.0   ;
  activeVar.resize(nbrIn+1) ; // update the active set
  activeVar[nbrIn] = varIn  ;
  nbrIn++                   ; // update the number of active variable
  isActiveVar[varIn] = 1    ;
}

void ACTIVE_SET::del_var(uword indVarOut) {
  betaActive.shed_row(indVarOut)        ; // update the vector of active parameters
  isActiveVar[activeVar[indVarOut]] = 0 ; // update the active set
  activeVar.shed_row(indVarOut)         ;
  nbrIn--                               ; // update the number of active variable
}

// must be sorted in increasing order
void ACTIVE_SET::add_vars(uword first_varIn, uword last_varIn) {
  for (int k=first_varIn; k<=last_varIn; k++) {
    add_var(k);
  }
}

void ACTIVE_SET::del_vars(uword first_indVarOut, uword last_indVarOut) {
  for (int k=last_indVarOut; k>=first_indVarOut; k--) {
    del_var(k);
    if (k == 0) {break;}
  }
}

// ______________________________________________________
// ACTIVE DATA CLASS
//
// Inherits from the ACTIVE_SET class
//
// Contains "activated" part of the data required for optimization

ACTIVE_DATA::ACTIVE_DATA(SEXP BETA0, REGRESSION_DATA data) :
  ACTIVE_SET(BETA0) {
  if (!betaActive.is_empty()) {
    xtxA = mat(data.xt * data.x.cols(activeVar)) ;
    smoothGrad = -data.xty + xtxA * betaActive   ;
    xAtxA = xtxA.rows(activeVar)                 ;
  }
}

void ACTIVE_DATA::add_var(uword varIn, REGRESSION_DATA data) {

  ACTIVE_SET::add_var(varIn) ; // update the active set

  vec new_col = data.xt * data.x.col(varIn);

  // if (lambda2 > 0) {
  // Adding the column corresponding to the structurating matrix
  // new_col += spS.col(var_in);
  // }

  // UPDATE THE xtxA AND xAtxA MATRICES
  if (!xAtxA.is_empty()) {
    xAtxA = join_cols(xAtxA, xtxA.row(varIn)) ;
  }
  xtxA  = join_rows(xtxA, new_col) ;
  xAtxA = join_rows(xAtxA, trans(xtxA.row(varIn))) ;
}

void ACTIVE_DATA::add_vars(uword first_varIn, uword last_varIn, REGRESSION_DATA data) {

  ACTIVE_SET::add_vars(first_varIn, last_varIn) ; // update the active set

  mat new_cols = data.xt * data.x.cols(first_varIn, last_varIn);

  // UPDATE THE xtxA AND xAtxA MATRICES
  if (!xAtxA.is_empty()) {
    xAtxA = join_cols(xAtxA, xtxA.rows(first_varIn, last_varIn)) ;
  }
  xtxA  = join_rows(xtxA, new_cols) ;
  xAtxA = join_rows(xAtxA, trans(xtxA.rows(first_varIn, last_varIn))) ;
}


void ACTIVE_DATA::del_var(uword indVarOut) {
  ACTIVE_SET::del_var(indVarOut) ;
  xtxA.shed_col(indVarOut)       ;
  xAtxA.shed_col(indVarOut)      ;
  xAtxA.shed_row(indVarOut)      ;
}

// indVarsOut MUST be sorted in the decreasing order
void ACTIVE_DATA::del_vars(uword first_indVarOut, uword last_indVarOut) {
  ACTIVE_SET::del_vars(first_indVarOut, last_indVarOut) ;
  xtxA.shed_cols(first_indVarOut, last_indVarOut)       ;
  xAtxA.shed_cols(first_indVarOut, last_indVarOut)      ;
  xAtxA.shed_rows(first_indVarOut, last_indVarOut)      ;
}

// handle vector of active variable, with method for adding/removing
// elements by group of variables


// GROUP_STRUCTURE::
// GROUP_STRUCTURE(SEXP PK):
//   pk (as<uvec>(PK)) // a vector with the group numbers
// {
//   K = pk.n_elem          ; // number of groups
//   gk_start.zeros(K)      ; // starting index for each group
//   gk_end.zeros  (K)      ; // starting index for each group
//   group.zeros(sum(pk))   ; // indicator of group belonging
//   is_active_grp.zeros(K) ;
//   nbrGrpIn = 0 ;

//   gk_start[0] = 0 ;
//   gk_end[0]   = pk[0];
//   for (int k = 1; k < K; k++) {
//     gk_start[k] = gk_start[k-1] + pk[k-1];
//     gk_end[k]   = gk_start[k] + pk[k];
//   }
//   for (int k = 0; k < K; k++) {
//     group.subvec(gk_start[k], gk_end[k]-1) = k*arma::ones<arma::uvec>(pk[k]);
//   }
// }

// void GROUP_STRUCTURE::add_var(uword grpIn) {

//   active_grp.resize(nbrGrpIn+1) ; // update the active set
//   active_grp[nbrGrpIn] = grpIn  ;
//   nbrGrpIn++                    ; // update the number of active variable
//   is_active_grp[grpIn] = 1      ;

//   for (int k=gk_start[grpIn]; k<gk_end[grpIn]; k++) {
//     ACTIVE_SET::add_var(k);
//   }
// }

// void GROUP_STRUCTURE::del_var(uword indGrpOut) {
//   is_active_grp[active_grp[indGrpOut]] = 0; // update the active set
//   active_grp.shed_row(indGrpOut)            ;
//   nbrGrpIn--                                ; // update the number of active variable

//   // Update set of active variables
//   uvec remove = sort(find(group.elem(activeVar) == active_grp[indGrpOut]),1) ;
//   for (int k=0; k<remove.n_elem; k++) {
//     ACTIVE_SET::del_var(remove[k]);
//   }
// }

