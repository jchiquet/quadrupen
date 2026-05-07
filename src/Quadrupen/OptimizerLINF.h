/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "Optimizer.h"
#include "PenaltyDense.h"
#include "ActiveSet.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

class OptimizerLINF : public Optimizer {

public:

  DensePenalty<DenseNorm::LINF> penalty_ ;
  OptimizerLINF() {} ; // needed
  OptimizerLINF(DensePenalty<DenseNorm::LINF>&, const List& control) ;

  uword quadratic_breg(
      vec& beta,
      vec &grad,
      const double& lambda,
      const vec& weights,
      RegressionData<mat> &data,
      uvec& unbounded,
      const double& accuracy,
      const uword& max_iter) ;

};


