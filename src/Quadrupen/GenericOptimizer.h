/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_GENERIC_OPTIMIZER_H
#define _quadrupen_GENERIC_OPTIMIZER_H

#include "RegressionData.h"
#include "ActiveSet.h"
#include "ActiveSetGroup.h"
#include "Penalty.h"

enum SolverType {FISTA, QUADRA};

#define ZERO 2e-16 // practical zero

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, Norm norm>
class GenericOptimizer {
public:

  GenericOptimizer() {} ;
  GenericOptimizer(Penalty<norm>&) ;

  uword conjugate_gradient(
      vec& x0,
      const mat& A,
      const vec& b,
      const double& accuracy,
      const uword& max_iter) ;
  
  uword fista(
      vec& beta,
      vec& grad,
      const double& lambda, 
      RegressionData<matrix> &data,
      ActiveSet<matrix>& set,
      const double& accuracy, 
      const uword& max_iter) ;

  uword fista_groupwise(
      vec& beta,
      const double& lambda, 
      RegressionData<matrix> &data,
      ActiveSetGroup<matrix>& set,
      const double& accuracy, 
      const uword& max_iter) ;
  
  // uword coordinate_descent(
  //     vec& beta,
  //     const double& lambda,
  //     ActiveSet<matrix>& set,
  //     const double& accuracy,
  //     const uword& max_iter) ;
      
protected:
  Penalty<norm> penalty_ ;

};

template <typename matrix, Norm norm>
inline GenericOptimizer<matrix, norm>::GenericOptimizer(
  Penalty<norm>& penalty) :
  penalty_ (penalty) {

}

template <typename matrix, Norm norm>
uword GenericOptimizer<matrix,norm>::conjugate_gradient(
    vec& x0,
    const mat& A,
    const vec& b,
    const double& accuracy,
    const uword& max_iter) {
  
  vec r = b - A * x0;
  vec p = r ;
  double rs_old = sum(square(r)) ;

  double rs_new = rs_old ;
  uword i = 0;
  double alpha ;
  mat Ap ;
  
  while ((sqrt(rs_new) > accuracy) & (i < max_iter)) {
    Ap = A * p;
    alpha = rs_old/dot(p,Ap) ;
    if (std::isfinite(alpha)) {
      x0 += alpha * p ;
      r -= alpha * Ap ;
      // Polak-Ribière update
      rs_new = dot(r,-alpha * Ap);
      p = r + rs_new/rs_old*p;
      rs_old = rs_new;
      i++;
    } else {
      break ;
    }
  }
  return(i) ;
}

template <typename matrix, Norm norm>
uword GenericOptimizer<matrix,norm>::fista(
    vec& beta,
    vec& grad,
    const double& lambda,
    RegressionData<matrix> &data,
    ActiveSet<matrix>& set,
    const double& accuracy,
    const uword& max_iter) {

  vec betak = beta  ; // output vector
  vec betal = beta  ;
  double delta = 2*accuracy  ; // change in beta
  double L = max( set.XATXA_.diag()) ; // Lipchitz constant

  double t0 = 1.0, tk ; // auxiliary variables in FISTA 
  uword iter = 0      ; // current iterate
  while ((delta > accuracy/beta.n_elem ) && (iter < max_iter)) {
    
    double l_num, l_den ;
    double f0, fk ;
    vec XATXA_betal = set.XATXA_ * betal ;
    f0 = dot(betal, .5 * XATXA_betal  - data.XTy_(set.A_)) ;
    grad((set.A_)) = - data.XTy_(set.A_) + XATXA_betal ;
    
    // Line search over L
    bool found=false;
    while(!found) {
      // Apply proximal operator (implemented in penalty object)
      vec prox_arg = betal - grad(set.A_)/L;
      betak = penalty_.proximal(prox_arg, lambda/L);

      fk = dot(betak, .5 * set.XATXA_ * betak - data.XTy_(set.A_)) ;
      l_num = 2 * (fk - f0 - dot(grad(set.A_), betak-betal));
      l_den = accu(pow(betak-betal,2));
      
      if ((L * l_den >= l_num) || (sqrt(l_den) < accuracy)) {
        found = true;
      } else {
        L = fmax(2*L, l_num/l_den);
      }

      R_CheckUserInterrupt();
    }

    // updating t
    tk = 0.5 * (1+sqrt(1+4*t0*t0));

    // updating s
    betal = betak + (t0-1)/tk * ( betak - beta );

    // preparing next iterate
    delta = sqrt(l_num);
    beta = betak;
    t0 = tk;
    iter++;

    R_CheckUserInterrupt();
  }

  return(iter) ;

}

template <typename matrix, Norm norm>
uword GenericOptimizer<matrix,norm>::fista_groupwise(
    vec& beta,
    const double& lambda,
    RegressionData<matrix> &data,
    ActiveSetGroup<matrix>& set,
    const double& accuracy,
    const uword& max_iter) {
  
  vec betak = beta  ; // output vector
  vec betal = beta  ;
  double delta = 2*accuracy  ; // change in beta
  double L = max( set.XATXA_.diag()) ; // Lipchitz constant
  
  double t0 = 1.0, tk ; // auxiliary variables in FISTA 
  uword iter = 0      ; // current iterate
  while ((delta > accuracy/beta.n_elem ) && (iter < max_iter)) {
    
    double l_num, l_den ;
    double f0, fk ;
    vec XATXA_betal = set.XATXA_ * betal ;
    f0 = dot(betal, .5 * XATXA_betal  - data.XTy_(set.A_)) ;
    vec grad = XATXA_betal - data.XTy_(set.A_);
    
    // Line search over L
    bool found=false;
    while(!found) {
      // Apply proximal operator (implemented in penalty object)
      vec prox_arg = betal - grad/L;
      betak = penalty_.proximal(prox_arg, lambda/L, set.grp_sizes_(set.G_));
      
      fk = dot(betak, .5 * set.XATXA_ * betak - data.XTy_(set.A_)) ;
      l_num = 2 * (fk - f0 - dot(grad, betak-betal));
      l_den = accu(pow(betak-betal,2));
      
      if ((L * l_den >= l_num) || (sqrt(l_den) < accuracy)) {
        found = true;
      } else {
        L = fmax(2*L, l_num/l_den);
      }
      
      R_CheckUserInterrupt();
    }
    
    // updating t
    tk = 0.5 * (1+sqrt(1+4*t0*t0));
    
    // updating s
    betal = betak + (t0-1)/tk * ( betak - beta );
    
    // preparing next iterate
    delta = sqrt(l_num);
    beta = betak;
    t0 = tk;
    iter++;
    
    R_CheckUserInterrupt();
  }
  
  return(iter) ;
  
}


// 
// template <typename matrix>
// uword Optimizer<matrix>::coordinate_descent(
//     vec& beta,
//     const double& lambda,
//     ActiveSet<matrix>& set,
//     mat& XTX,
//     const double& accuracy,
//     const uword& max_iter) {
//   

// int pathwise_enet(vec&  x0,
//                   mat& xtx,
//                   vec xty,
//                   vec& xtxw,
//                   double& pen,
//                   uvec &null,
//                   const double& gam   ,
//                   const double eps    ) {
  
//   double u, d               ; // temporary scalar
//   vec betak = beta         ; // output vector
//   double delta = 2*accuracy ; // change in beta
// 
//   double t0 = 1.0, tk ; // auxiliary variables in FISTA 
//   uword iter = 0      ; // current iterate
//   while ((delta > accuracy/beta.n_elem ) && (iter < max_iter)) {
// 
//     delta = 0;
//     for (uword j=0; j<beta.n_elem; j++) {
//       // Soft thresholding operator
//       u = beta(j) * (1+gam) + xty(j) - xtxw(j) ;
//       betak(j)  = fmax(1-pen/fabs(u),0) * u/(1+gam) ;
// 
//       // max(zeros(x.n_elem), 1-lambda/elt_norm(x)) % x;
//       
//       d = betak(j)-beta(j);
//       delta += pow(d,2);
//       xtxw  += d*xtx.col(j) ;
//     }
//     
//     // preparing next iterate
//     delta = sqrt(delta);
//     beta = betak;
//     iter++;
//     
//     R_CheckUserInterrupt();
//   }
//   
//   // null = sort(find(abs(betak) + (abs(-xty + xtxw) - pen) < ZERO), "descend") ;
//   return(iter);
// }


#endif

