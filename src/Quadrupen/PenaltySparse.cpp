/*
* Author: Julien CHIQUET
*         MIA Paris-Saclay
*/
  
#include "PenaltySparse.h"

using namespace Rcpp;
using namespace arma;

// ______________________________________________________
// L1 NORM A.K.A LASSO
template<>
vec SparsePenalty<SparseNorm::L1>::elt_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) % w);
}

template<>
vec SparsePenalty<SparseNorm::L1>::proximal(const vec& x, double lambda, const vec& w) {
  return(sign(x) % max(abs(x) - lambda * w, zeros<vec>(x.n_elem))) ;
}

// ______________________________________________________
// MCP

template<>
vec SparsePenalty<SparseNorm::MCP>::elt_norm(const vec& x, const vec& w, double lambda) {
  vec res = zeros<vec>(x.n_elem);
  for(uword i=0; i<x.n_elem; ++i) {
    double abs_xi = std::abs(x[i]);
    if (abs_xi <= eta_ * lambda) {
      res[i] = lambda * w[i] * abs_xi - (abs_xi * abs_xi) / (2.0 * eta_);
    } else {
      res[i] = 0.5 * eta_ * lambda * lambda * w[i]; // constant value (plateau)
    }
  }
  return res;
}

template<>
vec SparsePenalty<SparseNorm::MCP>::proximal(const vec& x, double lambda, const vec& w) {
  vec res = zeros<vec>(x.n_elem);
  for (uword i = 0; i < x.n_elem; ++i) {
    double abs_xi = std::abs(x[i]);
    double l = lambda * w[i];
    if (abs_xi <= l) {
      res[i] = 0.0;
    } else if (abs_xi <= eta_ * l) {
      res[i] = std::copysign(abs_xi - l, x[i]) / (1.0 - 1.0/eta_);
    } else {
      res[i] = x[i];
    }
  }
  return res;
}

// ______________________________________________________
// SCAD

template<>
vec SparsePenalty<SparseNorm::SCAD>::elt_norm(const vec& x, const vec& w, double lambda) {
  vec res = zeros<vec>(x.n_elem);
  for(uword i=0; i<x.n_elem; ++i) {
    double abs_xi = std::abs(x[i]);
    double l = lambda * w[i];
    if (abs_xi <= l) {
      res[i] = l * abs_xi;
    } else if (abs_xi <= eta_ * l) {
      res[i] = (2.0 * eta_ * l * abs_xi - abs_xi * abs_xi - l * l) / (2.0 * (eta_ - 1.0));
    } else {
      res[i] = (l * l * (eta_ + 1.0)) / 2.0;
    }
  }
  return res;
}

template<>
vec SparsePenalty<SparseNorm::SCAD>::proximal(const vec& x, double lambda, const vec& w) {
  vec res = zeros<vec>(x.n_elem);
  for (uword i = 0; i < x.n_elem; ++i) {
    double abs_xi = std::abs(x[i]);
    double l = lambda * w[i];
    if (abs_xi <= 2.0 * l) {
      // Soft-thresholding classique
      res[i] = std::copysign(std::max(0.0, abs_xi - l), x[i]);
    } else if (abs_xi <= eta_ * l) {
      // Zone intermédiaire (interpolation)
      res[i] = ((eta_ - 1.0) * x[i] - std::copysign(eta_ * l, x[i])) / (eta_ - 2.0);
    } else {
      res[i] = x[i]; // Pas de biais
    }
  }
  return res;
}

