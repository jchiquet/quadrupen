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

  Optimizer() ;
  Optimizer(RegressionData<matrix>&, Penalty&, SolverType&) ;
  
  int run(vec&, uvec&, vec&, mat&, double&, const double&, const uword&) ;
  
private:
  RegressionData<matrix> data_ ;
  Penalty penalty_ ;
  SolverType solver_ ;

  typedef int (Optimizer::*solver_ptr)(
      vec& beta,
      uvec& B,
      vec& grad,
      mat& XTX,
      double& lambda,
      const double& accuracy,
      const uword& max_iter) ;
  
  solver_ptr _current_solver_ptr ;

  int quadratic_breg(
      vec& beta,
      uvec& B,
      vec& grad,
      mat& XTX,
      double& lambda,
      const double& accuracy=1e-10,
      const uword& max_iter=50) ;
  
  int fista_breg(
      vec& beta,
      uvec& B,
      vec& grad,
      mat& XTX, 
      double& lambda, 
      const double& eps=1e-2, 
      const uword& max_iter=10000) ;
  
};

template <typename matrix>
Optimizer<matrix>::Optimizer(
  RegressionData<matrix>& data, Penalty& penalty,  SolverType& solver) :
  data_ (data), penalty_ (penalty), solver_(solver) {
  
  if (solver_ ==  QUADRA) {
    _current_solver_ptr  = &Optimizer::quadratic_breg ;
  } else if (solver == FISTA) {
    _current_solver_ptr  = &Optimizer::fista_breg ;
  }  
  
}

template <typename matrix>
int Optimizer<matrix>::run(
    vec& beta,
    uvec& B,
    vec& grad,
    mat& XTX,
    double& lambda,
    const double& accuracy,
    const uword& max_iter) 
{return((this->*_current_solver_ptr)(beta, B, grad, XTX, lambda, accuracy, max_iter));};

template <typename matrix>
int Optimizer<matrix>::quadratic_breg(
    vec& beta,
    uvec& B,
    vec& grad,
    mat& XTX,
    double& lambda,
    const double& accuracy,
    const uword& max_iter) {
  
  uvec all = regspace<uvec>(0,beta.n_elem-1) ;
  uvec toB = B ; // guys reaching the boundary after optimization
  
  grad = -data_.XTy_ + XTX * beta ;
  vec theta = -sign(grad.elem(B)) ;
  
  // SOLVE THE QUADRATIC PROBLEM
  uword iter = 0 ; // count the number of systems solved
  while ((toB.n_elem > 0) & (iter < max_iter)) {
    
    iter++;
    
    // Constructing the system (KKT)
    mat XX, XX_II ;
    vec XX_B = XTX.cols(B) * theta;
    uvec I = regspace<uvec>(0, data_.p_-1) ;
    I.shed_rows(B) ;
    
    if (I.is_empty()) {
      XX = sum(theta % XX_B.elem(B),0);
      double b  = (dot(theta, data_.XTy_.elem(B)) - lambda);
      beta.elem(B) = theta * (b/XX) ;
      // keep a trace of the current boundary
    } else {
      if (I.n_elem > 1) {
        XX_II = XTX.submat(I,I) ;
      } else {
        XX_II = XTX(I,I) ;
      }
      XX = join_rows(
        join_cols(sum(theta % XX_B.elem(B),0),XX_B.elem(I)),
        join_cols(strans(XX_B.elem(I)), XX_II));
      
      vec b = zeros<vec>(I.n_elem + 1) ;
      b[0] = dot(theta, data_.XTy_.elem(B)) - lambda;
      b.subvec(1,b.n_elem-1) = data_.XTy_.elem(I) ;
      
      // Solving via Cholesky factorization...
      mat R = chol(XX) ;
      vec tmp = solve(trimatu(R), solve(trimatl(strans(R)),b)) ;
      beta.elem(B) = theta * tmp[0] ;
      beta.elem(I) = tmp.subvec(1,tmp.n_elem-1) ;
    }
    // Handling guys reaching the boundary
    double bound = max(abs(beta.elem(B))); // current boundary
    toB = find(abs(beta) > bound);
    beta.elem(toB) = bound * sign(beta.elem(toB));
    B = unique(join_cols(B,toB));
    theta = sign(beta.elem(B)); // sign of the guys on the boundary
  }
  grad = -data_.XTy_ + XTX * beta ;
  
  // Handling guys leaving the boundary after optimization
  uvec toI  = find(abs(theta + sign(grad.elem(B))) > ZERO);
  B.shed_rows(toI) ;
  
  // If everyone is leaving the boundary, that's an issue...
  if (B.is_empty()) {
    throw std::runtime_error("Too much unstability");
  }
  
  return(iter) ;
}

template <typename matrix>
int Optimizer<matrix>::fista_breg(
    vec& beta,
    uvec& B,
    vec& grad,
    mat& XTX,
    double& lambda,
    const double& accuracy,
    const uword& max_iter) {
  
  colvec betak = beta  ; // output vector
  colvec betal = beta     ;
  uword iter = 0          ; // current iterate
  double delta = 2*accuracy  ; // change in beta
  double L0 = max( XTX.diag()) ; // Lipchitz constant
  
  double t0 = 1.0, tk;
  bool found=false;
  double f0, fk ;
  double l_num, l_den ;
  
  while ((delta > pow(accuracy,2)) && (iter < max_iter)) {
    
    f0 = as_scalar(.5 * strans(betal) * XTX * betal - strans(data_.XTy_) * betal) ;
    grad = -data_.XTy_ + XTX * betal ;
    
    // Line search over L
    while(!found) {
      betak = penalty_.proximal(betal - grad/L0, lambda /L0);
      
      fk = as_scalar(.5 * strans(betak) * XTX * betak - strans(data_.XTy_) * betak) ;
      l_num = as_scalar(2 * (fk - f0 - dot(grad, betak-betal) ));
      l_den = as_scalar(pow(norm(betak-betal,2),2));
      
      if ((L0 * l_den >= l_num) || (sqrt(l_den) < accuracy)) {
        found = true;
      } else {
        L0 = fmax(2*L0, l_num/l_den);
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
  
  // Evaluating the set of variable reaching the boundary
  B = find(abs(penalty_.elt_norm(beta) - penalty_.pen_norm(beta)) < ZERO );
  
  return(iter) ;
  
}

#endif

