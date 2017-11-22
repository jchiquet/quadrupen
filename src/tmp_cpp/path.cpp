#include "path.hpp"

using namespace Rcpp;
using namespace arma;

PATH::PATH (SEXP MAXFEAT) :
  maxFeat (as<uword>(MAXFEAT))
{}

void PATH::grid_penLevels(SEXP PENLEVELS,
			  SEXP PENLEN   ,
			  SEXP MINRATIO ,
			  double LMAX   ) {
  double lmax;
  if (PENLEVELS != R_NilValue) { // when the iser provide one, just import it
    penLevels = as<vec>(PENLEVELS)  ;
  } else { // otherwise compute a grid of penalties ona  log scale.
    penLevels = exp10(linspace(log10(LMAX), log10(as<double>(MINRATIO)*LMAX), as<uword>(PENLEN))) ;
  }
}
