/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _quadrupen_OPTIMIZER_H
#define _quadrupen_OPTIMIZER_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include "RegressionData_class.h"
#include "Penalty_class.h"

#define ZERO 2e-16 // practical zero

enum SolverType {FISTA, QUADRA};

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class Optimizer {
public:

  Optimizer() {} ;
  Optimizer(RegressionData<matrix>&, Penalty&, SolverType&) ;
  
  uword run(vec&, const double&, ActiveSet<matrix>&, mat&, const double&, const uword&) ;

  uword quadratic_breg(
      vec& beta,
      const double& lambda,
      ActiveSet<matrix>& set,
      mat& XTX,
      const double& accuracy=1e-10,
      const uword& max_iter=50) ;

  uword quadratic_enet(
      vec &beta0,
      const double &lambda ,
      ActiveSet<matrix> &set,
      mat&  XTX,
      const double& accuracy,
      const uword& max_iter) ;

  uword fista(
      vec& beta,
      const double& lambda, 
      ActiveSet<matrix>& set,
      mat& XTX, 
      const double& accuracy=1e-2, 
      const uword& max_iter=10000) ;

private:
  RegressionData<matrix> data_ ;
  Penalty penalty_ ;
  SolverType solver_ ;

  typedef uword (Optimizer::*solver_ptr)(
      vec& beta,
      const double& lambda,
      ActiveSet<matrix>& set,
      mat& XTX,
      const double& accuracy,
      const uword& max_iter) ;

  solver_ptr _current_solver_ptr ;

};

template <typename matrix>
Optimizer<matrix>::Optimizer(
  RegressionData<matrix>& data, Penalty& penalty,  SolverType& solver) :
  data_ (data), penalty_ (penalty), solver_(solver) {

  if (solver_ ==  QUADRA) {
    if (penalty_.type() == LINF) _current_solver_ptr = &Optimizer::quadratic_breg ;
    if (penalty_.type() == L1)   _current_solver_ptr = &Optimizer::quadratic_enet ;
  } else if (solver == FISTA) {
    _current_solver_ptr  = &Optimizer::fista ;
  }

}

template <typename matrix>
uword Optimizer<matrix>::run(
    vec& beta,
    const double& lambda,
    ActiveSet<matrix>& set,
    mat& XTX,
    const double& accuracy,
    const uword& max_iter)
{return((this->*_current_solver_ptr)(beta, lambda, set, XTX, accuracy, max_iter));};

template <typename matrix>
uword Optimizer<matrix>::quadratic_breg(
    vec& beta,
    const double& lambda,
    ActiveSet<matrix>& set,
    mat& XTX,
    const double& accuracy,
    const uword& max_iter) {

  uvec B = regspace<uvec>(0,beta.n_elem-1) ; B.shed_rows(set.A_) ;
  vec grad = -data_.XTy_ + XTX * beta ;
  vec theta = -sign(grad.elem(B)) ; // sign of the guys on the boundary

  uword iter = 0 ; // count the number of systems solved
  bool success = false ;
  while ((!success > 0) && (iter < max_iter)) {

    iter++;

    // SOLVE THE QUADRATIC PROBLEM
    vec XX_B = XTX.cols(B) * theta;
    if (set.A_.is_empty()) {
      double b  = (dot(theta, data_.XTy_.elem(B)) - lambda);
      beta.elem(B) = theta * (b/sum(theta % XX_B.elem(B),0)) ;
    } else {
      // Constructing the system (KKT)
      mat XX = join_rows(
        join_cols(sum(theta % XX_B.elem(B),0), XX_B.elem(set.A_)),
        join_cols(strans(XX_B.elem(set.A_)), set.XATXA_));
      vec b = join_cols(theta.t()*data_.XTy_.elem(B)-lambda, data_.XTy_.elem(set.A_)) ;
      // Solving via Cholesky factorization...
      mat R = chol(XX) ;
      vec tmp = solve(trimatu(R), solve(trimatl(strans(R)), b)) ;
      beta.elem(B) = theta * tmp[0] ;
      beta.elem(set.A_) = tmp.subvec(1,tmp.n_elem-1) ;
    }

    // Handling guys reaching the boundary (leaving the active set)
    double bound = penalty_.pen_norm(beta.elem(B)) ; // current boundary
    uvec ind_toB = find(abs(beta.elem(set.A_)) > bound) ;
    if (!ind_toB.is_empty()) {
      uvec toB = set.A_(ind_toB) ;
      set.del_vars(ind_toB) ;
      B = join_cols(B,toB);
      beta.elem(B) = bound * sign(beta.elem(B));
      theta = sign(beta.elem(B));
    } else {
      success = true ;
    }
  }

  // Guys leaving the boundary after optimization (activation)
  grad = -data_.XTy_ + XTX * beta ;
  uvec ind_toA  = find(theta == sign(grad.elem(B)));
  if (!ind_toA.is_empty()) {
    uvec toA = B(ind_toA) ;
    set.add_vars(toA, data_) ;
    B.shed_rows(ind_toA) ;
    if (B.is_empty()) {
      throw std::runtime_error("Every variable left the boudary. Try more regularization.");
    }
  }
  
  return(iter) ;
}

template <typename matrix>
uword Optimizer<matrix>::quadratic_enet(
    vec &beta0,
    const double &lambda ,
    ActiveSet<matrix> &set,
    mat&  XTX,
    const double& accuracy,
    const uword& max_iter) {

  uword iter = 1; // current iterate

  vec grad = - data_.XTy_(set.A_) + set.XATXA_ * beta0 ;
  vec theta = -sign(grad)   ; // vector of sign of the solution
  
  uvec A = find(abs(beta0) > ZERO) ; // vector of locally active variables
  theta.elem(A)   = sign(beta0.elem(A));

  // Solving the quadratic problem
  vec betak ;
  if (!set.R_.is_empty()) {
    betak = solve(trimatu(set.R_), 
                  solve(trimatl(strans(set.R_)), data_.XTy_(set.A_) - lambda * theta));
  } else {
    betak = solve(set.XATXA_, data_.XTy_(set.A_) - lambda*theta, 
                  solve_opts::likely_sympd + solve_opts::fast + 
                    solve_opts::force_approx) ;
  }

  // Check for swapping variables
  uvec swap = find(abs(sign(betak.elem(A)) - theta.elem(A)) > ZERO);
  uvec null ;
  if (swap.is_empty()) {
    null = swap ; // this is empty
    beta0 = betak;
  } else {
    swap = A.elem(swap);
    colvec beta0_swap = beta0.elem(swap);
    colvec betak_swap = betak.elem(swap);
    // first, go to zero for the swapped variable which cost the minimum
    vec eta = -beta0_swap / (betak_swap-beta0_swap);
    uword i_swap = eta.index_min();
    double scale = eta(i_swap) ;
    null = swap[i_swap];
    betak = beta0 + (betak-beta0) * scale ;
    // second, solve the problem after swaping the signs of the
    // incriminated variable
    beta0 = betak;
    beta0(null[0]) = -betak_swap[i_swap];

    A = find(abs(beta0) > ZERO) ; // new vector of active variables
    theta = -sign(grad)        ; // vector of sign of the solution
    theta.elem(A) = sign(beta0.elem(A));

    vec betal ;
    if (!set.R_.is_empty()) {
      betal = solve(trimatu(set.R_),
                    solve( trimatl(strans(set.R_)), data_.XTy_(set.A_) - lambda * theta));
    } else {
      betal = solve(set.XATXA_, data_.XTy_(set.A_) - lambda*theta, betak, 
                    solve_opts::likely_sympd + 
                      solve_opts::fast + 
                      solve_opts::force_approx) ;
    }
    iter++;
    
    // This is the gradient on the active part of the parameters
    grad = - data_.XTy_(set.A_) + set.XATXA_ * betal ;
    // if the sign is coherent, keep that one...
    if (fabs(grad(null[0]) + lambda * as_scalar(sign(betal(null)))) <= ZERO) {
      null = swap; // this is empty
      beta0 = betal ;
    } else {
      // otherwise, backtrack to betak
      beta0 = betak ;
    }
  }

  return(iter);
}


template <typename matrix>
uword Optimizer<matrix>::fista(
    vec& beta,
    const double& lambda,
    ActiveSet<matrix>& set,
    mat& XTX,
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
    f0 = as_scalar(.5 * strans(betal) * set.XATXA_ * betal - strans(data_.XTy_(set.A_)) * betal) ;
    vec grad = - data_.XTy_(set.A_) + set.XATXA_ * betal ;
    
    // Line search over L
    bool found=false;
    while(!found) {
      // Apply proximal operator (implemented in penalty object)
      betak = penalty_.proximal(betal - grad/L, lambda/L);

      fk = as_scalar(.5 * strans(betak) * set.XATXA_ * betak - strans(data_.XTy_(set.A_)) * betak) ;
      l_num = as_scalar(2 * (fk - f0 - dot(grad, betak-betal) ));
      l_den = as_scalar(pow(norm(betak-betal,2),2));

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
    found = false;
    iter++;

    R_CheckUserInterrupt();
  }

  return(iter) ;

}

#endif

