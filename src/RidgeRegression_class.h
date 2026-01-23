/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_RIDGE_H
#define _quadrupen_RIDGE_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "GenericRegularizer_class.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

class RidgeRegression : public GenericRegularizer{
 public:
   RidgeRegression(RegressionData&, const bool&, const List&);

  List solution_path(const mat&);
  
};

#endif

