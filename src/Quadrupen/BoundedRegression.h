/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _BoundedRegression_H
#define _BoundedRegression_H

#include "Regularizer.h"
#include "ActiveSet.h"
#include "PenaltySimple.h"
#include "OptimizerLINF.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

class BoundedRegression : public Regularizer<mat> {
public:

  // Specific to Bounded regression
  SimplePenalty<SimpleNorm::LINF> penalty_ ; // main penalty object 
  OptimizerLINF solver_ ; // Solvers for LINF penalty
// TODO: use a more simple forme of ActiveSet...
  ActiveSet<mat> set_   ; // Active set of variable and data
  vector<uvec> bounded_ ; // variables reaching the boundary (for all lambda values)
  
  BoundedRegression(RegressionData<mat>&, const List&, const List&);

  double get_lambda_max() {
    return(penalty_.dual_norm(data_.XTy_, lambda_factor_));
  }

  const sp_mat unbounded_var() { 
    vector<uword> rowA, colA ;
    uword current_col = 0;
    for (const auto& a : bounded_) {
      rowA.insert(rowA.end(), a.begin(), a.end()) ;
      colA.insert(colA.end(), a.n_elem, current_col) ;
      current_col++;
    }
    return sp_mat(join_cols(urowvec(rowA), urowvec(colA)),
                  vec(rowA.size(), fill::ones), data_.p_, bounded_.size(), true, false) ;
  }

  List solution_path(const List&);

  // Compute degrees of freedom for the current estimate
  double get_df() ; 
  
};

#endif

