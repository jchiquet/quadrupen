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

  
  // int quadra_enet(
  //     vec &beta,
  //     mat &R,
  //     mat &XTX,
  //     vec  xty,
  //     vec  sgn_grd,
  //     double &lambda,
  //     uvec   &null,
  //     const double& accuracy=1e-10) ;

  
};

template <typename matrix>
Optimizer<matrix>::Optimizer(
  RegressionData<matrix>& data, Penalty& penalty,  SolverType& solver) :
  data_ (data), penalty_ (penalty), solver_(solver) {

  if (solver_ ==  QUADRA) {
    if (penalty_.type() == LINF) _current_solver_ptr  = &Optimizer::quadratic_breg ;
    // if (penalty_.type() == L1) _current_solver_ptr  = &Optimizer::quadratic_enet ;
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

// template <typename matrix>
// int Optimizer<matrix>::quadra_enet(
//     vec &x0,
//     RegressionData& data,
//     ActiveSet &set,
//     vec  xty,
//     vec  sgn_grd,
//     double &lambda ,
//     uvec   &null,
//     double accuracy) {
//   
//   uword iter = 1; // current iterate
//   
//   uvec A = find(abs(x0) > ZERO) ; // vector of active variables
//   vec theta = -sgn_grd   ; // vector of sign of the solution
//   theta.elem(A)   = sign(x0.elem(A));
//   
//   // Solving the quadratic problem
//   vec x1 ;
//   if (usechol) {
//     x1 = solve(trimatu(R), solve(trimatl(strans(R)), xty - pen * theta));
//   } else {
//     x1 = solve(xAtxA, xty - pen*theta, solve_opts::fast + solve_opts::force_approx) ;
//   }
//   
//   // Check for swapping variables
//   uvec swap = find(abs(sign(x1.elem(A)) - theta.elem(A)) > ZERO);
//   if (swap.is_empty()) {
//     null = swap ; // this is empty
//     x0 = x1;
//   } else {
//     swap = A.elem(swap);
//     colvec x0_swap = x0.elem(swap);
//     colvec x1_swap = x1.elem(swap);
//     // first, go to zero for the swapped variable which cost the minimum
//     vec gamma = -x0_swap / (x1_swap-x0_swap);
//     uword i_swap = gamma.index_min();
//     double scale = gamma(i_swap) ;
//     null = swap[i_swap];
//     x1 = x0 + (x1-x0) * scale ;
//     // second, solve the problem after swaping the signs of the
//     // incriminated variable
//     x0 = x1;
//     x0(null[0]) = -x1_swap[i_swap];
//     
//     A = find(abs(x0) > ZERO) ; // vector of active variables
//     theta = -sgn_grd        ; // vector of sign of the solution
//     theta.elem(A)   = sign(x0.elem(A));
//     
//     vec x2 ;
//     if (usechol) {
//       x2 = solve(trimatu(R), solve( trimatl(strans(R)), xty - pen * theta));
//     } else {
//       x2 = cg(xAtxA, xty - pen*theta, x1, tol) ;
//     }
//     iter++;
//     
//     // This is the gradient on the active part of the parameters
//     vec grd = -xty + xAtxA * x2;
//     // if the sign is coherent, keep that one...
//     if (fabs(grd(null[0]) + pen * as_scalar(sign(x2(null)))) <= ZERO) {
//       null = swap; // this is empty
//       x0 = x2 ;
//     } else {
//       // otherwise, backtrack to x1
//       x0 = x1 ;
//     }
//   }
//   
//   return(iter);
// }


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
  double L = max( XTX.diag()) ; // Lipchitz constant

  double t0 = 1.0, tk ; // auxiliary variables in FISTA 
  uword iter = 0      ; // current iterate
  while ((delta > accuracy/beta.n_elem ) && (iter < max_iter)) {
    
    double l_num, l_den ;
    double f0, fk ;
    f0 = as_scalar(.5 * strans(betal) * XTX * betal - strans(data_.XTy_) * betal) ;
    vec grad = - data_.XTy_ + XTX * betal ;
    
    // Line search over L
    bool found=false;
    while(!found) {
      // Apply proximal operator (implemented in penalty object)
      betak = penalty_.proximal(betal - grad/L, lambda/L);

      fk = as_scalar(.5 * strans(betak) * XTX * betak - strans(data_.XTy_) * betak) ;
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

  // Evaluating the set of active variables (between the boundary)
  // A = find(abs(penalty_.elt_norm(beta) - penalty_.pen_norm(beta)) >= ZERO );

  return(iter) ;

}

#endif

