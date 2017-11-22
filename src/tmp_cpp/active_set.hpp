/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_ACTIVE_SET_H
#define _quadrupen_ACTIVE_SET_H

#include <RcppArmadillo.h>
#include "data_reg.hpp"

using namespace Rcpp;
using namespace arma;

class ACTIVE_SET {

protected:
  // VARIABLES FOR HANDLING THE ACTIVE SET
  vec  betaActive  ; // vector of currently activated features
  uvec activeVar   ; // set of currently activated variables
  uvec isActiveVar ; // a vector to check if a variable is already in the active set
  uword nbrIn      ; // number of active variables

public:
  ACTIVE_SET (SEXP BETA0);

  // ACTIVE SET HANDLING
  void add_var(uword varIn)      ; // add a variable in the active set
  void del_var(uword indVarOut)  ; // remove the variable activated in postition ind_var_out
  void add_vars(uword first_varIn    , uword last_varIn    ) ; // the same for a set 
  void del_vars(uword first_indVarOut, uword last_indVarOut) ; // of contiguous variables

  const vec & get_betaActive() const { return betaActive ; }
  const uvec & get_activeVar() const { return activeVar  ; }
  bool is_activeVar(uword i) { return (isActiveVar(i) == 1) ; }
};


class ACTIVE_DATA: public ACTIVE_SET {
private:
  bool l2Reg      ; // TODO
  mat xAtxA       ;
  mat xtxA        ;
  vec smoothGrad  ;

public:
  ACTIVE_DATA(SEXP BETA0, REGRESSION_DATA data) ;
  void add_var (uword varIn, REGRESSION_DATA data) ; // add a variable
  void del_var(uword indVarOut)                    ; // remove a variable
  void add_vars(uword first_varIn, uword last_varIn, REGRESSION_DATA data) ; // the same for a set 
  void del_vars(uword first_indVarOut, uword last_indVarOut) ; // of contiguous variables

  const mat& GramMat() const { return xAtxA ; }
} ;

#endif


// class ACTIVE_SET_GROUPWISE: public ACTIVE_SET {

// protected:
//   // VARIABLES FOR HANDLING THE ACTIVE SET
//   uvec active_grp    ; // set of currently activated groups
//   uvec is_active_grp ; // a vector to check if a variable is already in the active set
//   uword nbrGrpIn     ; // number of active variables
//   uvec  pk           ; // vector of group numbers
//   uword K            ; // number of groups
//   uvec  gk_start     ; // starting index for each group
//   uvec  gk_end       ; // starting index for each group
//   uvec  group        ; // indicator of group belonging for each variable

// public:
//   ACTIVE_SET_GROUPWISE (SEXP BETA0, SEXP PK);

//   // ACTIVE SET HANDLING
//   void add_var(uword grpIn)   ; // add a group of variables in the active set
//   void del_var(uword indGrpOut) ; // remove the group of variable activated in postition ind_var_class

//   const uvec   & get_active_grp() const { return active_grp ; }

// };

// // ACTIVE_ALGORITHM
// // vec GROUP_LASSO::get_dist2opt(double lambda) {

// //   vec normGrad = max(zeros(K), grp_norm(smooth_grad) - lambda) ;
// //   // normGrad.elem(active_group) = grp_norm(smooth_grad + lambda * betaA);

// //   return(normGrad) ;
// // }

