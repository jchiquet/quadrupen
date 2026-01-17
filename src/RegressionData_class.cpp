/*
 * Author: Julien CHIQUET, INRAE
 *         julien.chiquet@inrae.fr
 *         Statistique et Génome
 *         MIA Paris-Saclay
 */

#include "RegressionData_class.h"

using namespace Rcpp;
using namespace arma;

RegressionData::RegressionData(
  const Environment& dataModel,
  const bool& center, const bool& scale):
  n_       (as<uword>  (dataModel["n"])) , // sample size
  p_       (as<uword>  (dataModel["d"])) , // # of features
  X_       (as<mat>    (dataModel["X"])) , // matrix of features
  y_       (as<vec>    (dataModel["y"])) , // response vector
  S_       (as<sp_mat> (dataModel["S"])) , // structuring matrix
  weights_ (as<vec>    (dataModel["wy"]))  // observation weights
{
  standardize(center, scale) ;
}

void RegressionData::standardize(const bool center, const bool scale) {

  // TODO: OBSERVATION WEIGHTS

  if (center) {
    X_bar_ = mean(X_, 0).t();
    y_bar_ = mean(y_) ;
    // do not center X to preserve the sparsity structure
  } else {
    X_bar_ = zeros(p_) ;
    y_bar_ = 0;
  }

  if (scale) {
    norm_X_ = sqrt(trans(sum(square(X_))) - n_ * square(X_bar_));
    for (uword i=0; i<p_; i++) {
      X_.col(i) /= norm_X_(i);
    }
    X_bar_ /= norm_X_ ;
  } else {
    norm_X_ = ones(p_);
  }
  norm_y_ = sqrt(sum(square(y_))) ;

  if (center) {
    XTy_ = X_.t() * (y_- y_bar_) ;
    for (uword i=0; i<p_; i++) {
      XTy_(i) -=  sum(y_-y_bar_) * X_bar_(i);
    }
  } else {
    XTy_ = X_.t() * y_ ;
  }
  
  centered_ = center ;
  scaled_   = scale  ;

}

