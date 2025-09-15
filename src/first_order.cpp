/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#include "first_order.h"

int fista_lasso(vec  &x0   ,
		mat  &xtx  ,
		vec xty  ,
		double &pen,
		uvec   &null,
		double &L0 ,
		double eps) {

  colvec xk = x0  ; // output vector
  colvec  s = x0  ;
  int iter = 0, j = 0   ; // current iterate
  double delta = 2*eps  ; // change in beta
  int max_iter  = 10000 ; // max. number of iteration
  double L  = L0        ; // Lipschitz's constant

  double t0 = 1.0, tk;
  bool found=false;
  double f0, fk ;
  colvec df ;

  double l_num, l_den ;

  while ((delta > eps/x0.n_elem ) && (iter < max_iter)) {

    f0 = as_scalar(.5 * strans(s) * xtx * s - strans(xty) * s) ;
    df = - xty + xtx * s ;

    // Line search over L
    while(!found) {
      // apply proximal operator of the Lasso
      xk = s - df/L ;
      for (j=0; j<x0.n_elem; j++) {
        xk(j) = fmax(0,1-(pen/L)/fabs(xk(j))) * xk(j);
      }

      fk = as_scalar(.5 * strans(xk) * xtx * xk - strans(xty) * xk) ;
      l_num = as_scalar(2 * (fk - f0 - dot(df, xk-s) ));
      l_den = as_scalar(pow(norm(xk-s,2),2));

      if ((L * l_den >= l_num) || (sqrt(l_den) < eps)) {
	      found = true;
      } else {
	L = fmax(2*L, l_num/l_den);
      }

      R_CheckUserInterrupt();
    }

    // updating t
    tk = 0.5 * (1+sqrt(1+4*t0*t0));

    // updating s
    s = xk + (t0-1)/tk * ( xk - x0 );

    // preparing next iterate
    delta = sqrt(l_num);
    x0 = xk;
    t0 = tk;
    found = false;
    iter++;

    R_CheckUserInterrupt();
  }

  null = sort(find(abs(xk) + (abs(df) - pen) < ZERO), "descend") ;

  return(iter);
}

int fista_breg(vec    &x0,
	       const mat &xtx,
	       const vec &xty,
	       vec    &grd,
	       double &pen,
	       double &L0 ,
	       double eps) {

  colvec xk = x0        ; // output vector
  colvec  s = x0        ;
  int iter = 0          ; // current iterate
  double delta = 2*eps  ; // change in beta
  int max_iter  = 10000 ; // max. number of iteration

  double t0 = 1.0, tk;
  bool found=false;
  double f0, fk ;
  colvec lbd (x0.n_elem) ;
  lbd.fill(pen);

  double l_num, l_den ;

  while ((delta > eps*eps) && (iter < max_iter)) {

    f0 = as_scalar(.5 * strans(s) * xtx * s - strans(xty) * s) ;
    grd = - xty + xtx * s ;

    // Line search over L
    while(!found) {
      xk = proximal_inf(s - grd/L0, lbd /L0);

      fk = as_scalar(.5 * strans(xk) * xtx * xk - strans(xty) * xk) ;
      l_num = as_scalar(2 * (fk - f0 - dot(grd, xk-s) ));
      l_den = as_scalar(pow(norm(xk-s,2),2));

      if ((L0 * l_den >= l_num) || (sqrt(l_den) < eps)) {
	found = true;
      } else {
	L0 = fmax(2*L0, l_num/l_den);
      }

      R_CheckUserInterrupt();
    }

    // updating t
    tk = 0.5 * (1+sqrt(1+4*t0*t0));

    // updating s
    s = xk + (t0-1)/tk * ( xk - x0 );

    // preparing next iterate
    delta = sqrt(l_num);
    x0 = xk;
    t0 = tk;
    found = false;
    iter++;

    R_CheckUserInterrupt();
  }

  return(iter) ;
}

int fista_grp(vec   &x0,
	      uvec  &pk,
	      mat   &xtx,
	      vec   xty,
	      vec   &pen,
	      uword &grpnorm,
	      double &L0    ,
	      double eps   ) {

  colvec xk = x0  ; // output vector
  colvec  s = x0  ;
  int iter = 0, j = 0      ; // current iterate
  double delta = 2*eps  ; // change in beta
  int max_iter  = 10000 ; // max. number of iteration
  double L  = L0        ; // Lipschitz's constant

  double t0 = 1.0, tk;
  bool found=false;
  double f0, fk ;
  colvec df ;

  double l_num, l_den ;

  while ( (delta > eps/x0.n_elem) && (iter < max_iter)) {

    f0 = as_scalar(.5 * strans(s) * xtx * s - strans(xty) * s) ;
    df = (xtx * s - xty) ;

    // Line search over L
    while(!found) {
      // apply proximal operator of the group-Lasso
      if (grpnorm == -1) { // -1 is our code for "inf" norm
	xk = proximal_grpinf(s - df/L, pk, pen, L);
      } else {
	xk = proximal_grp2(s - df/L, pk, pen, L);
      }

      fk = as_scalar(.5 * strans(xk) * xtx * xk - strans(xty) * xk) ;
      l_num = as_scalar(2 * (fk - f0 - dot(df, xk-s) ));
      l_den = as_scalar(pow(norm(xk-s,2),2));

      if (L * l_den >= l_num  | sqrt(l_den) < eps) {
	found = true;
      } else {
	L = fmax(2*L, l_num/l_den);
      }

      R_CheckUserInterrupt();
    }

    // updating t
    tk = 0.5 * (1+sqrt(1+4*t0*t0));

    // updating s
    s = xk + (t0-1)/tk * ( xk - x0 );

    // preparing next iterate
    delta = sqrt(l_num);
    x0 = xk;
    t0 = tk;
    found = false;
    iter++;

    R_CheckUserInterrupt();
  }

  return(iter);
}

vec proximal_grp2(vec u,
		  uvec pk,
		  vec  lambda,
		  double L) {

  int j,k,ind = 0 ;

  vec grp_norm2 = grp_norm(u, pk, 2, 1);
  for (k=0; k<pk.n_elem; k++) {
    for (j=ind; j<(ind+pk.at(k)); j++) {
      u(j) = fmax(0,1-(lambda(j)/L)/grp_norm2(j)) * u(j);
    }
    ind += pk.at(k);
  }

  return(u);
}

vec proximal_grpinf(vec u,
		    uvec pk,
		    vec lambda,
		    double L) {

  int ind = 0 ;
  for (int k=0; k<pk.n_elem; k++) {
    u.subvec(ind,ind+pk.at(k)-1) = proximal_inf(u.subvec(ind,ind+pk.at(k)-1),lambda.subvec(ind,ind+pk.at(k)-1)/L);
    ind += pk.at(k);
  }
  return(u);
}


vec proximal_inf(vec v,
		 vec lambda) {

  int p = v.n_elem;
  vec u, proj;
  vec out = zeros<vec>(p);

  if ( as_scalar(sum(abs(v) / lambda)) >= 1) {

    // Reordonnons les valeurs absolues
    u = sort(abs(v), "descend");

    // valeurs des coordonnées projetees si non nulles (problème dual)
    proj = (cumsum(u) - lambda)/linspace<vec>(1,p,p);

    // selection des coordonnees non nulles (problème dual)
    uvec maxs = sort(find(u-proj>ZERO), "descend") ;
    double thresh = proj[maxs[0]];

    // solution du programme primal
    // On garde les valeurs les plus petites, et on seuille le reste à une valeur commune +- thresh
    for (int k=0; k < p;k++) {
      if (fabs(v(k)) > ZERO) {
	if (v(k) > 0) {
	  out(k) = fmin(fabs(v(k)),thresh);
	} else {
	  out(k) = -fmin(fabs(v(k)),thresh);
	}
      }
    }
  }
  return(out);
}

int pathwise_enet(vec&  x0,
		  mat& xtx,
		  vec xty,
		  vec& xtxw,
		  double& pen,
		  uvec &null,
		  double& gam   ,
		  double eps    ) {

  colvec xk = x0 ; // output vector
  int j, iter  = 0     ; // current iterate
  int max_iter = 10000 ; // max. number of iteration
  double delta = 2*eps ; // change in beta
  double u, d          ; // temporary scalar

  while ((delta > eps/x0.n_elem ) && (iter < max_iter)) {

    delta = 0;
    for (j=0; j<x0.n_elem; j++) {
      // Soft thresholding operator
      u = x0(j) * (1+gam) + xty(j) - xtxw(j) ;
      xk(j)  = fmax(1-pen/fabs(u),0) * u/(1+gam) ;
      d = xk(j)-x0(j);
      delta += pow(d,2);
      xtxw  += d*xtx.col(j) ;
    }

    // preparing next iterate
    delta = sqrt(delta);
    x0 = xk;
    iter++;

    R_CheckUserInterrupt();
  }

  null = sort(find(abs(xk) + (abs(-xty + xtxw) - pen) < ZERO), "descend") ;
  return(iter);
}

// int pathwise_grpl1linf(vec  &x0,
// 		       uvec &pk,
// 		       mat &xtx,
// 		       vec xty,
// 		       vec &xtxw,
// 		       vec &pen,
// 		       double &gam   ,
// 		       double eps    ) {

//   colvec xk = x0 ; // output vector
//   int j, k, iter = 0     ; // current iterate
//   int ind      = 0     ; // current variable
//   int max_iter = 10 ; // max. number of iteration
//   double delta = 2*eps ; // change in beta
//   double normu       ;
//   vec u          ;
//   vec d          ;
//   vec us          ;

//   while ((delta > eps/x0.n_elem ) && (iter < max_iter)) {

//     delta = 0;
//     ind   = 0;
//     // blockwise coordinate descent
//     for (k=0; k<pk.n_elem; k++) {

//       u = x0.subvec(ind,ind+pk.at(k)-1) * (1+gam) + xty.subvec(ind,ind+pk.at(k)-1) - xtxw.subvec(ind,ind+pk.at(k)-1) ;

//       us = sort(abs(u),1);
//       if (as_scalar(sum(abs(u))) > pen(ind)) {

// 	// valeurs des coordonnees projetees si non nulles
// 	vec proj = (cumsum(us) - pen(ind)) /linspace<vec>(1,x0.n_elem,x0.n_elem);

// 	// selection des coordonnees non nulles
// 	uword jmax;
// 	proj.max(jmax);

// 	for (j=ind; j<(ind+pk.at(k)); j++) {

// 	  if (abs(u(iter)) > eps/100) {
// 	    if (u(iter) > 0) {
// 	      xk(j) = fmin(abs(u(i)),proj(jmax));
// 	    } else {
// 	      xk(j) = -fmin(abs(u(i)),proj(jmax));
// 	    }
// 	  }
// 	}
//       } else {
// 	xk.subvec(ind,ind+pk.at(k)-1) = zeros<vec>(pk.at(k));
//       }
//       d = xk.subvec(ind,ind+pk.at(k)-1)-x0.subvec(ind,ind+pk.at(k)-1);
//       delta += norm(d,2);
//       xtxw  += xtx.cols(ind,ind+pk.at(k)-1) * d;

//       ind += pk.at(k);

//     }

//     // preparing next iterate
//     delta = sqrt(delta);
//     x0 = xk;
//     iter++;

//     R_CheckUserInterrupt();
//   }

//   return(iter);
// }
