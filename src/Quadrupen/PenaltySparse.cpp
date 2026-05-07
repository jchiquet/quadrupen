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
vec SparsePenalty<SparseNorm::L1>::elt_dual_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) / w);
}

template<>
double SparsePenalty<SparseNorm::L1>::pen_norm(const vec& x, const vec& w, double lambda) {
  return(accu(elt_norm(x, w)));
}

template<>
double SparsePenalty<SparseNorm::L1>::dual_norm(const vec& x, const vec& w, double lambda) {
  return(max(elt_dual_norm(x, w))) ;
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
    if (abs_xi <= gamma_ * lambda) {
      res[i] = lambda * w[i] * abs_xi - (abs_xi * abs_xi) / (2.0 * gamma_);
    } else {
      res[i] = 0.5 * gamma_ * lambda * lambda * w[i]; // constant value (plateau)
    }
  }
  return res;
}

template<>
double SparsePenalty<SparseNorm::MCP>::dual_norm(const vec& x, const vec& w, double lambda) {
  return arma::max(arma::abs(x) / w);
}

template<>
vec SparsePenalty<SparseNorm::MCP>::proximal(const vec& x, double lambda, const vec& w) {
  vec res = zeros<vec>(x.n_elem);
  for (uword i = 0; i < x.n_elem; ++i) {
    double abs_xi = std::abs(x[i]);
    double l = lambda * w[i];
    if (abs_xi <= l) {
      res[i] = 0.0;
    } else if (abs_xi <= gamma_ * l) {
      res[i] = std::copysign(abs_xi - l, x[i]) / (1.0 - 1.0/gamma_);
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
    } else if (abs_xi <= gamma_ * l) {
      res[i] = (2.0 * gamma_ * l * abs_xi - abs_xi * abs_xi - l * l) / (2.0 * (gamma_ - 1.0));
    } else {
      res[i] = (l * l * (gamma_ + 1.0)) / 2.0;
    }
  }
  return res;
}

// Comme MCP, la condition d'entrée (dual) est dominée par la pente à l'origine
// Donc comme le Lasso
template<>
double SparsePenalty<SparseNorm::SCAD>::dual_norm(const vec& x, const vec& w, double lambda) {
  return arma::max(arma::abs(x) / w);
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
    } else if (abs_xi <= gamma_ * l) {
      // Zone intermédiaire (interpolation)
      res[i] = ((gamma_ - 1.0) * x[i] - std::copysign(gamma_ * l, x[i])) / (gamma_ - 2.0);
    } else {
      res[i] = x[i]; // Pas de biais
    }
  }
  return res;
}

