/*
* Author: Julien CHIQUET
*         MIA Paris-Saclay
*/
  
#include "PenaltySimple.h"

using namespace Rcpp;
using namespace arma;

// ______________________________________________________
// L1 NORM A.K.A LASSO
template<>
vec SimplePenalty<SimpleNorm::L1>::elt_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) % w);
}

template<>
vec SimplePenalty<SimpleNorm::L1>::elt_dual_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) / w);
}

template<>
double SimplePenalty<SimpleNorm::L1>::pen_norm(const vec& x, const vec& w, double lambda) {
  return(accu(elt_norm(x, w)));
}

template<>
double SimplePenalty<SimpleNorm::L1>::dual_norm(const vec& x, const vec& w, double lambda) {
  return(max(elt_dual_norm(x, w))) ;
}

template<>
vec SimplePenalty<SimpleNorm::L1>::proximal(const vec& x, double lambda, const vec& w) {
  return(sign(x) % max(abs(x) - lambda * w, zeros<vec>(x.n_elem))) ;
}

// ______________________________________________________
// LINF NORM A.K.A BOUNDED REGRESSION
template<>
vec SimplePenalty<SimpleNorm::LINF>::elt_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) % w);
}

template<>
vec SimplePenalty<SimpleNorm::LINF>::elt_dual_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) / w);
}

template<>
double SimplePenalty<SimpleNorm::LINF>::pen_norm(const vec& x, const vec& w, double lambda) {
  return(max(elt_norm(x, w)));
}

template<>
double SimplePenalty<SimpleNorm::LINF>::dual_norm(const vec& x, const vec& w, double lambda) {
  return(sum(elt_dual_norm(x,w))) ;
}

template<>
vec SimplePenalty<SimpleNorm::LINF>::proximal(const vec& x, double lambda, const vec& w) {
  uword p = x.n_elem;
  vec abs_x_w = arma::abs(x) / w;

  if (accu(abs_x_w) <= lambda) {
    return zeros<vec>(p);
  }

  // Project onto the l1 ball

  // Reordering absolute values
  vec u = sort(abs_x_w, "descend");
  
  // values of the projected coordinate if non zero (dual problem)
  // vec proj = (cumsum(u) - lambda)/linspace<vec>(1,p,p);
  // 
  // find critical index
  // uword i = max(find(u - proj >= 0)) ;

  double cumsum_u_w = 0.0;
  double cumsum_w_sq = 0.0;
  double rho = 0.0;
  
  for (uword j = 0; j < p; ++j) {
    
    cumsum_u_w += u(j) * (w(j) * w(j)); //
    cumsum_w_sq += w(j) * w(j);
    
    double t = (cumsum_u_w - lambda) / cumsum_w_sq;
    
    if (j < p - 1 && u(j + 1) <= t) {
      rho = t;
      break;
    }
    if (j == p - 1) rho = t;
  }
  
  vec res = sign(x) % min(abs(x), rho * w );

  return(res);
}

// ______________________________________________________
// MCP

template<>
vec SimplePenalty<SimpleNorm::MCP>::elt_norm(const vec& x, const vec& w, double lambda) {
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
double SimplePenalty<SimpleNorm::MCP>::dual_norm(const vec& x, const vec& w, double lambda) {
  return arma::max(arma::abs(x) / w);
}

template<>
vec SimplePenalty<SimpleNorm::MCP>::proximal(const vec& x, double lambda, const vec& w) {
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
vec SimplePenalty<SimpleNorm::SCAD>::elt_norm(const vec& x, const vec& w, double lambda) {
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
double SimplePenalty<SimpleNorm::SCAD>::dual_norm(const vec& x, const vec& w, double lambda) {
  return arma::max(arma::abs(x) / w);
}

template<>
vec SimplePenalty<SimpleNorm::SCAD>::proximal(const vec& x, double lambda, const vec& w) {
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



// ______________________________________________________
// L2 NORM SQUARED A.K.A RIDGE
template<>
vec SimplePenalty<SimpleNorm::L2>::elt_norm(const vec& x, const vec& w, double lambda) {
  return(arma::pow(x, 2) % w);
}

template<>
vec SimplePenalty<SimpleNorm::L2>::elt_dual_norm(const vec& x, const vec& w, double lambda) {
  return(arma::pow(x, 2) / w);
}

template<>
double SimplePenalty<SimpleNorm::L2>::pen_norm(const vec& x, const vec& w, double lambda) {
  return(sum(elt_norm(x, w)));
}

// Slight abuse: take the Fenchel conjugate of L2^2
template<>
double SimplePenalty<SimpleNorm::L2>::dual_norm(const vec& x, const vec& w, double lambda) {
  return(.25*sum(elt_norm(x, w))) ;
}

template<>
vec SimplePenalty<SimpleNorm::L2>::proximal(const vec& x, double lambda, const vec& w) {
  return(x / (1+2*lambda*w));
}
