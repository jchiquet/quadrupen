/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_SIMPLE_OPTIMIZER_H
#define _quadrupen_SIMPLE_OPTIMIZER_H

#include "Optimizer.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix, SimpleNorm norm> class SimpleOptimizer: public Optimizer<matrix> {
public:
  
  using Optimizer<matrix>::algorithm_  ;
  using Optimizer<matrix>::accuracy_   ;
  using Optimizer<matrix>::maxiter_    ;
  using Optimizer<matrix>::maxfeat_    ;
  using Optimizer<matrix>::verbosity_  ;
  using Optimizer<matrix>::monitoring_ ;
  using Optimizer<matrix>::iter_       ;
  using Optimizer<matrix>::inner_iter_ ;
  using Optimizer<matrix>::gap_        ;
  using Optimizer<matrix>::J_          ;
  using Optimizer<matrix>::D_          ;
  using Optimizer<matrix>::J_vec_      ;
  using Optimizer<matrix>::D_vec_      ;

  SimpleOptimizer() {} ;
  SimpleOptimizer(SimplePenalty<norm>&, const List&) ;

  SimplePenalty<norm> penalty_ ;
  
  void optimality_gap(
      vec& beta,
      vec& grad,
      const double& lambda,
      const double& gamma, 
      RegressionData<matrix> &data,
      ActiveSet<matrix>& set, uword type) ;
  
  uword solve(
      vec& beta,
      vec& grad,
      const double& lambda, 
      const double& gamma, 
      RegressionData<matrix> &data,
      ActiveSet<matrix>& set
  ) ;
  
  uword fista(
      vec& beta,
      vec& grad,
      const double& lambda, 
      RegressionData<matrix> &data,
      ActiveSet<matrix>& set,
      const double& accuracy, 
      const uword& max_iter
  ) ;
  
  virtual uword quadratic(
      vec &beta,
      vec &grad,
      const double &lambda ,
      RegressionData<matrix> &data,
      ActiveSet<matrix> &set,
      const double& accuracy,
      const uword& max_iter) = 0 ;
  
};

template <typename matrix, SimpleNorm norm>
SimpleOptimizer<matrix, norm>::SimpleOptimizer(
    SimplePenalty<norm>& penalty, const List& control) : 
  penalty_ (penalty), Optimizer<matrix>(control) {}

template <typename matrix, SimpleNorm norm>
uword SimpleOptimizer<matrix,norm>::fista(
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

template <typename matrix, SimpleNorm norm>
void SimpleOptimizer<matrix,norm>::optimality_gap(
    vec& beta,
    vec& grad,
    const double& lambda,
    const double& gamma,
    RegressionData<matrix> &data,
    ActiveSet<matrix>& set, uword type) {
  
  // nu equals the max |gradient|
  double nu = arma::norm(grad, "inf");
  double loss = .5 * pow(data.norm_y_ ,2) + 
    dot(beta, .5 * set.XATXA_ * beta - data.XTy_(set.A_)) ;
  double old_J = J_, old_D = D_ ;
  J_ = loss - dot(beta, grad(set.A_))  ;
  uvec Ac ;
  
  switch (type) {
  case 1: // Grandvalet's bound
    Ac = find(grad > nu); // set of adversarial variables outside the boundary
    D_ = J_ * (1 - lambda/nu) - 
      (pow(lambda,2)/(2*gamma))*((lambda*(data.p_-Ac.n_elem))/nu + 
      pow(arma::norm(grad(Ac),2)/nu,2)-data.p_);
    break;
  case 2: // Fenchel's bound
    if (nu < lambda) nu = lambda;
    D_ = loss * (1+pow(lambda/nu,2)) + sum(abs(lambda*beta)) + 
      (lambda/nu)*(dot(beta,data.XTy_(set.A_))-pow(data.norm_y_,2));
    break;
  default: 
    D_ = datum::inf ;
  break;
  }
  
  // keep the smallest bound reached so far for a given lambda value
  if ((old_J < J_) && (old_D - D_) < (old_J - J_)) {D_ = old_D ; }
  
}

template <typename matrix, SimpleNorm norm>
uword SimpleOptimizer<matrix,norm>::solve(
    vec& beta,
    vec& grad,
    const double& lambda,
    const double& gamma, 
    RegressionData<matrix> &data,
    ActiveSet<matrix>& set) {
  
  // variable associated with the highest violation of KKT conditions 
  vec optimality = penalty_.elt_norm(grad) - lambda ;
  uword status = 0, var_in = optimality.index_max() ;
  gap_ = std::max(0.0, optimality(var_in)) ;
  
  iter_ = 0 ; bool success = true ; 
  J_ = datum::inf ; D_ = datum::inf ;
  while ((gap_ > accuracy_) && (iter_ <= maxiter_)) {
    R_CheckUserInterrupt();
    iter_++;
    
    // VARIABLE ACTIVATION IF APPLICABLE
    if (set.is_in_[var_in] == 0) { // Is var_in already in the active set?
      set.add_var(var_in, data) ;
      beta.resize(beta.size()+1) ; // update the vector of active parameters
      beta.tail(1) = - 1e-3 * sign(grad(var_in)) ;
      if (verbosity_) {Rprintf("\tnewly added variable %i\n",var_in);}
    } else if (verbosity_) {Rprintf("\talready in %i\n",var_in);}
    
    // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    if (algorithm_ == FISTA) {
      inner_iter_.push_back(
        fista(beta, grad, lambda, data, set, 1e-10, 10000)
      );
    } else { // QUADRA solver
      try {
        inner_iter_.push_back(
          quadratic(beta, grad, lambda, data, set, 1e-5, 10000)
        );
      } catch (std::runtime_error& error) {
        if (verbosity_ > 0) {
          Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
        }
        success = false ;
      }
    }
    
    // OPTIMALITY TESTING
    grad = - data.XTy_ + set.XTXA_ * beta ;
    optimality = penalty_.elt_norm(grad) - lambda ;
    var_in = optimality.index_max() ;
    gap_ = std::max(0.0, optimality(var_in)) ;
    
    if (monitoring_ > 0) {
      optimality_gap(beta, grad, lambda, gamma, data, set, monitoring_) ;
      J_vec_.push_back(J_) ;
      D_vec_.push_back(D_) ;
    }
    
  }
  if (verbosity_) Rprintf("\tcurrent gap = %f\n",gap_) ;
  
  // Checking convergence status
  if (iter_ >= maxiter_)     { status = 1 ; }
  if (set.size() > maxfeat_) { status = 2 ; }
  if (!success)              { status = 3 ; }
  
  return status ;
}

#endif

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

