/*
 * Author: Julien CHIQUET
 */
#include "quadratic.h"

int quadra_enet(vec &x0,
		mat &R,
		mat &xAtxA,
		vec  xty,
		vec  sgn_grd,
		double &pen ,
		uvec   &null,
		bool   usechol,
		double tol) {

  uword iter = 1; // current iterate

  uvec A = find(abs(x0) > ZERO) ; // vector of active variables
  vec theta = -sgn_grd   ; // vector of sign of the solution
  theta.elem(A)   = sign(x0.elem(A));

  // Solving the quadratic problem
  vec x1 ;
  if (usechol) {
    x1 = solve(trimatu(R), solve(trimatl(strans(R)), xty - pen * theta));
  } else {
    x1 = cg(xAtxA, xty - pen*theta, x0, tol) ;
  }

  // Check for swapping variables
  uvec swap = find(abs(sign(x1.elem(A)) - theta.elem(A)) > ZERO);
  if (swap.is_empty()) {
    null = swap ; // this is empty
    x0 = x1;
  } else {
    swap = A.elem(swap);
    colvec x0_swap = x0.elem(swap);
    colvec x1_swap = x1.elem(swap);
    // first, go to zero for the swapped variable which cost the minimum
    vec gamma = -x0_swap / (x1_swap-x0_swap);
    uword i_swap = gamma.index_min();
    double scale = gamma(i_swap) ;
    null = swap[i_swap];
    x1 = x0 + (x1-x0) * scale ;
    // second, solve the problem after swaping the signs of the
    // incriminated variable
    x0 = x1;
    x0(null[0]) = -x1_swap[i_swap];

    A = find(abs(x0) > ZERO) ; // vector of active variables
    theta = -sgn_grd        ; // vector of sign of the solution
    theta.elem(A)   = sign(x0.elem(A));

    vec x2 ;
    if (usechol) {
      x2 = solve(trimatu(R), solve( trimatl(strans(R)), xty - pen * theta));
    } else {
      x2 = cg(xAtxA, xty - pen*theta, x1, tol) ;
    }
    iter++;

    // This is the gradient on the active part of the parameters
    vec grd = -xty + xAtxA * x2;
    // if the sign is coherent, keep that one...
    if (fabs(grd(null[0]) + pen * as_scalar(sign(x2(null)))) <= ZERO) {
      null = swap; // this is empty
      x0 = x2 ;
    } else {
      // otherwise, backtrack to x1
      x0 = x1 ;
    }
  }
  
  return(iter);
}
