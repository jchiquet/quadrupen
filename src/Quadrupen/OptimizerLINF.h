/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_OPTIMIZER_LINF_H
#define _quadrupen_OPTIMIZER_LINF_H

#include "Optimizer.h"
#include "PenaltySimple.h"
#include "ActiveSet.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

class OptimizerLINF : public Optimizer {

public:

  SimplePenalty<SimpleNorm::LINF> penalty_ ;
  OptimizerLINF() {} ;
  OptimizerLINF(SimplePenalty<SimpleNorm::LINF>&, const List& control) ;

  uword quadratic_breg(
      vec& beta,
      vec &grad,
      const double& lambda,
      const vec& weights,
      RegressionData<mat> &data,
      uvec& active,
      const double& accuracy,
      const uword& max_iter) ;

  mat updateCholeskyFromExisting(const mat& R, const vec& b) ;
  
};

#endif

