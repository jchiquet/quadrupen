/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "GenericRegularizer_class.h"

using namespace Rcpp;
using namespace arma;

GenericRegularizer::~GenericRegularizer() {} ;

GenericRegularizer::GenericRegularizer(
  RegressionData& data, const bool& intercept, const List& regParam) :
  data_ (data), intercept_ (intercept)
{
  penalty_ = Penalty() ;
}

void GenericRegularizer::lambda_seq(const List& regParam) {
  if (regParam[0] != R_NilValue) {
    lambda_  = as<vector<double>>(regParam["lambda"]) ;
  } else {
    double lambda_max = this->get_lambda_max() ;
    lambda_ = conv_to<vector<double>>::from(
      logspace(
        log10(lambda_max),
        log10(as<double>(regParam["min_ratio"])*lambda_max),
        as<uword>(regParam["n_lambda"])
      )
    );
  }
}
