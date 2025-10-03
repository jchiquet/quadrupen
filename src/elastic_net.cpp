/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "quadrupen_headers.hpp"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List elastic_net_cpp(
     SEXP BETA0                 , //
		 SEXP X                     , // regressor matrix
		 const arma::vec y          , // response vector
		 const arma::sp_mat Struct  , // Structuring matrix
		 SEXP LAMBDA1               , // 
		 double n_lambda            , //
		 const double min_ratio     , //
		 const arma::vec penscale   , // penalty weights
		 const double lambda2       , // the smooth (ridge) penalty
		 const bool intercept       , // boolean for intercept mode
		 const bool normalize       , // boolean for standardizing the predictor
		 const arma::vec weights    , // observation weights (not use at the moment)
		 const bool naive           , // naive elastic-net or not
		 const double eps           , // precision required
		 const arma::uword max_iter , // max # of iterates of the active set
		 const arma::uword max_feat , // max # of variables activated
		 const arma::uword fun      , // solver (0=quadra, 1=pathwise, 2=fista)
		 const arma::uword verbose  , // int for verbose mode (0/1/2)
		 const bool sparse          , // boolean for sparse mode
		 const bool usechol         , // use Cholesky decomposition or not
		 const arma::uword monitor    // convergence monitoring (1 == Grandvalet's bound ;-) 2 == Fenchel duality gap)
		 ) {

  vec    xty   ; // responses to predictors vector
  vec    xbar  ; // mean of the predictors
  vec    normx ; // norm of the predictors
  double normy ; // norm of the response
  double ybar  ; // mean of the response
  const double eps2 = pow(eps, 2) ;
  const uword n(y.n_elem)      ; // sample size
  const uword p(Struct.n_cols)      ; // problem size

  mat x        ;
  mat xt       ;
  sp_mat sp_x  ;
  sp_mat sp_xt ;
  if (sparse == 1) { // Check how x is encoded for reading
    sp_x = as<sp_mat>(X) ;
    standardize(sp_x, y, intercept, normalize, penscale, xty, normx, normy, xbar, ybar) ;
    sp_xt = sp_x.t() ;
  } else {
    x = as<mat>(X) ;
    standardize(x, y, intercept, normalize, penscale, xty, normx, normy, xbar, ybar) ;
    xt = x.t();
  }

  // STRUCTURATING MATRIX
  sp_mat diag_S = spdiags(sqrt(lambda2)*pow(penscale,-1/2), ivec({0}), p, p) ;
  const sp_mat S = diag_S * Struct * diag_S  ; // sparsely encoded structuring matrix
  mat SAA ; // densely encoded active counterpart

  // VECTOR OF TUNING PARAMETER FOR THE L1-PENALTY
  vec lambda1 = get_lambda1(LAMBDA1, n_lambda, min_ratio, max(abs(xty)));
  n_lambda = lambda1.n_elem  ; // # of penalty levels

  // Initializing "first level" variables (outside of the lambda1 loop)
  mat  R                                 ; // Cholesky decomposition of XAtXA
  uvec A                                 ; // set of currently activated variables
  vec  betaA                             ; // vector of currently activated parameters
  mat  xtxA                              ; // t(x) * x_A  covariance matrix
  mat  xAtxA                             ; // t(x_A) * x_A + S covariance matrix of the activated variable plus SAA matrix
  vec  xtxw                              ; // t(x_A) * x_A * beta(A)
  vec  grd       = -xty                  ; // smooth part of the gradient
  vec  mu        = zeros<vec>(n_lambda)  ; // the intercept term
  vec  max_grd   = zeros<vec>(n_lambda)  ; // a vector with the successively reach duality gap
  vec  converge  = zeros<vec>(n_lambda)  ; // a vector indicating if convergence occured (0/1/2)
  uvec it_active = zeros<uvec>(n_lambda) ; // # of loop in the active set for each lambda1
  uvec it_optim                          ; // # of loop in the optimization process for each loop of the active se
  double L0      = 1.0 + lambda2         ; // Lipschitz constant for proximal methods
  vec  timing      (n_lambda)            ; // succesive timing in
  vec  df          (n_lambda)            ; // degrees of freedom
  wall_clock timer                       ; // clock

  // Initializing "second level" variables (within the active set - for a fixed value of lamdba)
  uword var_in                           ; // currently added variable
  uword nbr_in   = 0                     ; // # of currently added variables
  uword nbr_opt  = 0                     ; // # of current calls to the optimization routine
  uvec  are_in   = zeros<uvec>(p)        ; // a vector to check if a variable is already in the active set
  List  out_optim                        ; // the list of output of the optimization function
  bool  success_optim = true             ; // was the internal system resolution successful?
  uvec  null                             ; // stores the variables which go to zero during optimization
  vec   grd_norm (p)                     ; // current value of the grd_norm for each variable
  mat   nonzeros                         ; // contains non-zero value of beta
  mat   iA                               ; // contains row indices of the non-zero values
  mat   jA                               ; // contains column indices of the non-zero values

  vec beta0 ;
  // WARM START
  if (BETA0 != R_NilValue) {
    beta0 = as<vec>(BETA0) ;
    A = find(beta0 != 0) ;
    betaA = beta0.elem(A) ;
    if (sparse == 1) {
      // WRONG - DO IT THE RIGHT WAY
      xtxA = mat(sp_xt * sp_x.col(0)) ;
    } else {
      xtxA = mat(xt * x.cols(A)) ;
    }
    if (lambda2 > 0) {
      for (uword i=0; i<A.n_elem;i++) {
	      xtxA.col(i) = xtxA.col(i) + S.col(A(i));
	      are_in(A(i)) = 1;
      }
    }
    grd += xtxA * betaA    ;
    nbr_in = A.n_elem      ;
    xAtxA = xtxA.rows(A)   ;
    if ((fun == 0) & (usechol)) {
      R = chol(xAtxA) ;
    }
    if (fun == 1) {
      xtxw(nbr_in) = dot(xAtxA.col(nbr_in),betaA);
    }
  }

  // Additional variable for convergence monitoring
  vec D_hat    ;
  vec D_star   ;
  vec J_hat    ;
  mat J_star   ;

  // _____________________________________________________________
  //
  // START THE LOOP OVER LAMBDA
  timer.tic();
  for (int m=0; m<n_lambda; m++) {
    if (verbose == 2) {
      Rprintf("\n lambda1 = %f",lambda1(m)) ;
      Rprintf("\n nb active variables = %i",nbr_in) ;
    }
    // _____________________________________________________________
    //
    // START THE ACTIVE SET ALGORITHM
    // _____________________________________________________________
    //

    // dual norm of gradient for unactive variable
    grd_norm = abs(grd) - lambda1[m] ;
    // gradient for active variables
    grd_norm.elem(A) = abs(grd.elem(A) + lambda1[m] * sign(betaA)) ;
    // variable associated with the highest optimality violation
    var_in = grd_norm.index_max() ;

    max_grd[m] = grd_norm(var_in) ;
    if (max_grd[m] < 0) {max_grd[m] = 0;}

    while ((max_grd[m] > eps) && (it_active[m] < max_iter)) {
      // _____________________________________________________________
      //
      // (1) VARIABLE ACTIVATION IF APPLICABLE
      // _____________________________________________________________

      // Check if the variable is already in the active set
      if (are_in[var_in] == 0) {
	      if (sparse == 1) {
	        add_var_enet(n, nbr_in, var_in, betaA, A, sp_x, sp_xt, xtxA, xAtxA, xtxw, R, lambda2, xbar, S, usechol, fun) ;
	      } else {
	        add_var_enet(n, nbr_in, var_in, betaA, A, x, xt, xtxA, xAtxA, xtxw, R, lambda2, xbar, S, usechol, fun) ;
	      }
	      if (verbose == 2) {Rprintf("newly added variable %i\n",var_in);}
	      are_in[var_in] = 1;
	      nbr_in++;
      } else {
	      if (verbose == 2) {Rprintf("already in %i\n",var_in);}
      }

      // _____________________________________________________________
      //
      // (2) OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      // _____________________________________________________________
      //

      it_optim.reshape(nbr_opt + 1,1) ;
      switch (fun) {
        case 1 :
	        it_optim[nbr_opt] = pathwise_enet(betaA, xAtxA, xty.elem(A), xtxw, lambda1[m], null, lambda2, eps2);
	        break;
        case 2 :
	        it_optim[nbr_opt] = fista_lasso(betaA, xAtxA, xty.elem(A), lambda1[m], null, L0, eps2);
	        break;
        default:
	        try {
	          it_optim[nbr_opt] = quadra_enet(betaA, R, xAtxA, xty.elem(A), sign(grd.elem(A)), lambda1[m], null, usechol, eps);
	        } catch (std::runtime_error &error) {
	          if (verbose > 0) {
	          Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
	          }
	          success_optim = false ;
	      }
      }
      // update the smooth part of the gradient
      grd = -xty + xtxA * betaA;
      nbr_opt++;

      // _____________________________________________________________
      //
      // (3) VARIABLE DELETION IF APPLICABLE
      // _____________________________________________________________
      //
      // removing variables zeroed during optimization
      if (!null.is_empty()) {
	      if (verbose == 2) {
	        for (uword j=0; j<null.n_elem; j++) {Rprintf("removing variable %i\n",null[j]);}
	      }
	      remove_var_enet(nbr_in,are_in,betaA,A,xtxA,xAtxA,xtxw,R,null,usechol,fun) ;
      }

      // _____________________________________________________________
      //
      // (4) OPTIMALITY TESTING
      // _____________________________________________________________

      // dual norm of gradient for unactive variable
      grd_norm = abs(grd) - lambda1[m] ;
      // dual norm of gradient for active variables
      grd_norm.elem(A) = abs(grd.elem(A) + lambda1[m] * sign(betaA)) ;
      // variable associated with the highest optimality violation
      var_in  = grd_norm.index_max() ;
      max_grd[m]  = grd_norm(var_in) ;
      if (max_grd[m] < 0) {max_grd[m] = 0;}

      if (monitor > 0) {
	      // _____________________________________________________________
	      //
	      // (OPTIONAL) FOLLOWING CONVERGENCE BY COMPLETE MONITORING
	      // _____________________________________________________________
	      bound_to_optimal(betaA, xAtxA, xty, grd, lambda1[m], lambda2, normy, A, monitor, J_hat, D_hat) ;
      }

      // Moving to the next iterate
      it_active[m]++;

      R_CheckUserInterrupt();
    }

    // degrees of freedom
    df[m] = get_df_enet(lambda2, R, xAtxA, S, A, fun);

    // the reference parameter (obtained once optimum is met)
    if (monitor > 0) {
      if (it_active[m] > 0) {
	      J_star = join_cols(J_star, ones(it_active[m],1) * J_hat[nbr_opt-1]) ;
      }
    }

    // Record the time ellapsed
    timing[m] = timer.toc() ;

    // Checking convergence status
    if (it_active[m] >= max_iter) {
      converge[m] = 1;
    }
    if (nbr_in > max_feat) {
      converge[m] = 2 ;
    }
    if (!success_optim) {
      converge[m] = 3;
    }

    // Stop now if relevant
    if (converge[m] == 2 || converge[m] == 3) {
      lambda1     =    lambda1.subvec(0,m-1) ;
      converge    =  converge.subvec(0,m)    ;
      max_grd     =   max_grd.subvec(0,m-1)  ;
      it_active   = it_active.subvec(0,m)    ;
      timing      =    timing.subvec(0,m)    ;
      df          =    df.subvec(0,m)        ;
      break;
    } else {
      nonzeros = join_cols(nonzeros, betaA/(normx.elem(A) % penscale.elem(A)));
      iA = join_cols(iA, m*ones(betaA.n_elem,1) );
      jA = join_cols(jA, conv_to<mat>::from(A) ) ;
      if (intercept == 1) {
	      mu[m] = dot(betaA, xbar.elem(A)) ;
      }
    }
  }
  // END OF THE LOOP OVER LAMBDA
  if (!naive) {
    nonzeros *= 1+lambda2;
    mu = ybar - (1+lambda2) * mu;
  } else {
    mu = ybar - mu;
  }

  // Updating monitored quantities
  if (monitor > 0) {
    D_star = J_hat - J_star;
  }

  return List::create(Named("nzeros")     = nonzeros ,
		      Named("iA")         = iA       ,
		      Named("jA")         = jA       ,
		      Named("mu")         = mu       ,
		      Named("normx")      = normx    ,
		      Named("lambda1")    = lambda1  ,
		      Named("df")         = df       ,
		      Named("nbr.in")     = nbr_in   ,
		      Named("it.active")  = it_active,
		      Named("it.optim")   = it_optim ,
		      Named("max.grd")    = max_grd  ,
		      Named("timing")     = timing   ,
		      Named("delta.hat")  = D_hat    ,
		      Named("delta.star") = D_star   ,
		      Named("converge")   = converge );

}


// =========================================================================== //
//
// 
// V2
// [[Rcpp::export]]
Rcpp::List elastic_net2_cpp(
    SEXP BETA0                   , // initial vector of coefficients
    const Environment &dataModel , // data structure
       arma::vec lambda1         , // vector of L1 penalties
    const double lambda2         , // scalar for the amount L2 penalty
    const List control             // config of the optimisation 
) {

  const uword n             = dataModel["n"]  ; // sample size
  const uword p             = dataModel["d"]  ; // problem size
  const SEXP &X             = dataModel["X"]  ; // design matrix
  const arma::vec &y        = dataModel["y"]  ; // response vector
  const arma::sp_mat& S     = dataModel["S"]  ; // Structuring matrix
  const arma::vec &penscale = dataModel["wx"] ;  // penalty weights
  const bool &intercept     = dataModel["has_intercept"] ; // boolean for intercept mode
  const arma::vec &xty      = dataModel["xty"]    ; // responses to predictors vector
  const arma::vec xbar      = dataModel["mean_X"] ; // mean of the predictors
  const arma::vec &normx    = dataModel["norm_X"] ; // norm of the predictors
  const double ybar         = dataModel["mean_y"] ; // mean of the predictors
  const double normy        = dataModel["norm_y"] ; // norm of the response
  const arma::vec& weights  = dataModel["wy"]     ; // observation weights (not use at the moment)
  const bool sparse         = dataModel["sparse_encoding"] ; // boolean for sparse mode
  
  const double eps           = control["threshold"] ; // precision required
  const arma::uword max_iter = control["max.iter"]  ; // max # of iterates of the active set
  const arma::uword max_feat = control["max.feat"]  ; // max # of variables activated
  const arma::uword fun      = control["method"]    ; // solver (0=quadra, 1=pathwise, 2=fista)
  const arma::uword verbose  = control["verbose"]   ; // int for verbose mode (0/1/2)
  const bool usechol         = control["usechol"]   ; // use Cholesky decomposition or not
  const bool naive           = control["naive"]     ; // use Cholesky decomposition or not
  const arma::uword monitor  = control["monitor"]   ; // convergence monitoring (1 == Grandvalet's bound ;-) 2 == Fenchel duality gap)

  const double eps2 = pow(eps, 2) ;
  
  mat x        ;
  mat xt       ;
  sp_mat sp_x  ;
  sp_mat sp_xt ;
  if (sparse) { // Check how x is encoded for reading
    sp_x = as<sp_mat>(X) ;
    sp_xt = sp_x.t() ;
  } else {
    x = as<mat>(X) ;
    xt = x.t();
  }

  // Initializing "first level" variables (outside of the lambda1 loop)
  uword n_lambda = lambda1.n_elem        ; // # of penalty levels
  mat  R                                 ; // Cholesky decomposition of XAtXA
  uvec A                                 ; // set of currently activated variables
  vec  betaA                             ; // vector of currently activated parameters
  mat  xtxA                              ; // t(x) * x_A  covariance matrix
  mat  xAtxA                             ; // t(x_A) * x_A + S covariance matrix of the activated variable plus SAA matrix
  vec  xtxw                              ; // t(x_A) * x_A * beta(A)
  vec  grd       = -xty                  ; // smooth part of the gradient
  vec  mu        = zeros<vec>(n_lambda)  ; // the intercept term
  vec  max_grd   = zeros<vec>(n_lambda)  ; // a vector with the successively reach duality gap
  vec  converge  = zeros<vec>(n_lambda)  ; // a vector indicating if convergence occurred (0/1/2)
  uvec it_active = zeros<uvec>(n_lambda) ; // # of loop in the active set for each lambda1
  uvec it_optim                          ; // # of loop in the optimization process for each loop of the active se
  double L0      = 1.0 + lambda2         ; // Lipschitz constant for proximal methods
  vec  timing      (n_lambda)            ; // successive timing in
  vec  df          (n_lambda)            ; // degrees of freedom
  wall_clock timer                       ; // clock
  
  // Initializing "second level" variables (within the active set - for a fixed value of lambda1)
  uword var_in                           ; // currently added variable
  uword nbr_in   = 0                     ; // # of currently added variables
  uword nbr_opt  = 0                     ; // # of current calls to the optimization routine
  uvec  are_in   = zeros<uvec>(p)        ; // a vector to check if a variable is already in the active set
  List  out_optim                        ; // the list of output of the optimization function
  bool  success_optim = true             ; // was the internal system resolution successful?
  uvec  null                             ; // stores the variables which go to zero during optimization
  vec   grd_norm (p)                     ; // current value of the grd_norm for each variable
  mat   nonzeros                         ; // contains non-zero value of beta
  mat   iA                               ; // contains row indices of the non-zero values
  mat   jA                               ; // contains column indices of the non-zero values
  
  vec beta0 ;
  // WARM START
  if (BETA0 != R_NilValue) {
    beta0 = as<vec>(BETA0) ;
    A = find(beta0 != 0) ;
    betaA = beta0.elem(A) ;
    if (sparse) {
      // WRONG - DO IT THE RIGHT WAY
      xtxA = mat(sp_xt * sp_x.col(0)) ;
    } else {
      xtxA = mat(xt * x.cols(A)) ;
    }
    if (lambda2 > 0) {
      for (uword i=0; i<A.n_elem;i++) {
        xtxA.col(i) = xtxA.col(i) + S.col(A(i));
        are_in(A(i)) = 1;
      }
    }
    grd += xtxA * betaA    ;
    nbr_in = A.n_elem      ;
    xAtxA = xtxA.rows(A)   ;
    if ((fun == 0) & (usechol)) {
      R = chol(xAtxA) ;
    }
    if (fun == 1) {
      xtxw(nbr_in) = dot(xAtxA.col(nbr_in),betaA);
    }
  }
  
  // Additional variable for convergence monitoring
  vec D_hat, D_star, J_hat ; mat J_star ;
  
  // _____________________________________________________________
  //
  // START THE LOOP OVER LAMBDA
  timer.tic();
  for (uword m=0; m<n_lambda; m++) {
    if (verbose == 2) {
      Rprintf("\n lambda1 = %f", lambda1(m)) ;
      Rprintf("\n nb active variables = %i",nbr_in) ;
    }
    // _____________________________________________________________
    //
    // START THE ACTIVE SET ALGORITHM
    // _____________________________________________________________
    //
    
    // dual norm of gradient for inactive variables
    grd_norm = abs(grd) - lambda1[m] ;
    // gradient for active variables
    grd_norm.elem(A) = abs(grd.elem(A) + lambda1[m] * sign(betaA)) ;
    // variable associated with the highest violation of optimality conditions 
    var_in = grd_norm.index_max() ;
    
    max_grd[m] = grd_norm(var_in) ;
    if (max_grd[m] < 0) {max_grd[m] = 0;}
    
    while ((max_grd[m] > eps) && (it_active[m] < max_iter)) {
      // _____________________________________________________________
      //
      // (1) VARIABLE ACTIVATION IF APPLICABLE
      // _____________________________________________________________
      
      // Check if the variable is already in the active set
      if (are_in[var_in] == 0) {
        if (sparse) {
          add_var_enet(n, nbr_in, var_in, betaA, A, sp_x, sp_xt, xtxA, xAtxA, xtxw, R, lambda2, xbar, S, usechol, fun) ;
        } else {
          add_var_enet(n, nbr_in, var_in, betaA, A, x, xt, xtxA, xAtxA, xtxw, R, lambda2, xbar, S, usechol, fun) ;
        }
        if (verbose == 2) {Rprintf("newly added variable %i\n",var_in);}
        are_in[var_in] = 1;
        nbr_in++;
      } else {
        if (verbose == 2) {Rprintf("already in %i\n",var_in);}
      }
      
      // _____________________________________________________________
      //
      // (2) OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      // _____________________________________________________________
      //
      
      it_optim.reshape(nbr_opt + 1,1) ;
      switch (fun) {
      case 1 :
        it_optim[nbr_opt] = pathwise_enet(betaA, xAtxA, xty.elem(A), xtxw, lambda1[m], null, lambda2, eps2);
        break;
      case 2 :
        it_optim[nbr_opt] = fista_lasso(betaA, xAtxA, xty.elem(A), lambda1[m], null, L0, eps2);
        break;
      default:
        try {
          it_optim[nbr_opt] = quadra_enet(betaA, R, xAtxA, xty.elem(A), sign(grd.elem(A)), lambda1[m], null, usechol, eps);
        } catch (std::runtime_error &error) {
          if (verbose > 0) {
            Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
          }
          success_optim = false ;
        }
      }
      // update the smooth part of the gradient
      grd = -xty + xtxA * betaA;
      nbr_opt++;
      
      // _____________________________________________________________
      //
      // (3) VARIABLE DELETION IF APPLICABLE
      // _____________________________________________________________
      //
      // removing variables zeroed during optimization
      if (!null.is_empty()) {
        if (verbose == 2) {
          for (uword j=0; j<null.n_elem; j++) {Rprintf("removing variable %i\n",null[j]);}
        }
        remove_var_enet(nbr_in,are_in,betaA,A,xtxA,xAtxA,xtxw,R,null,usechol,fun) ;
      }
      
      // _____________________________________________________________
      //
      // (4) OPTIMALITY TESTING
      // _____________________________________________________________
      
      // dual norm of gradient for unactive variable
      grd_norm = abs(grd) - lambda1[m] ;
      // dual norm of gradient for active variables
      grd_norm.elem(A) = abs(grd.elem(A) + lambda1[m] * sign(betaA)) ;
      // variable associated with the highest optimality violation
      var_in  = grd_norm.index_max() ;
      max_grd[m]  = grd_norm(var_in) ;
      if (max_grd[m] < 0) {max_grd[m] = 0;}
      
      if (monitor > 0) {
        // _____________________________________________________________
        //
        // (OPTIONAL) FOLLOWING CONVERGENCE BY COMPLETE MONITORING
        // _____________________________________________________________
        bound_to_optimal(betaA, xAtxA, xty, grd, lambda1[m], lambda2, normy, A, monitor, J_hat, D_hat) ;
      }
      
      // Moving to the next iterate
      it_active[m]++;
      
      R_CheckUserInterrupt();
    }
    
    // degrees of freedom
    df[m] = get_df_enet(lambda2, R, xAtxA, S, A, fun);
    
    // the reference parameter (obtained once optimum is met)
    if (monitor > 0) {
      if (it_active[m] > 0) {
        J_star = join_cols(J_star, ones(it_active[m],1) * J_hat[nbr_opt-1]) ;
      }
    }
    
    // Record the time ellapsed
    timing[m] = timer.toc() ;
    
    // Checking convergence status
    if (it_active[m] >= max_iter) converge[m] = 1 ;
    if (nbr_in > max_feat)        converge[m] = 2 ;
    if (!success_optim)           converge[m] = 3 ;

    // Stop now if relevant
    if (converge[m] == 2 || converge[m] == 3) {
      lambda1     =    lambda1.subvec(0,m-1) ;
      converge    =  converge.subvec(0,m)    ;
      max_grd     =   max_grd.subvec(0,m-1)  ;
      it_active   = it_active.subvec(0,m)    ;
      timing      =    timing.subvec(0,m)    ;
      df          =    df.subvec(0,m)        ;
      break;
    } else {
      nonzeros = join_cols(nonzeros, betaA/(normx.elem(A) % penscale.elem(A)));
      iA = join_cols(iA, m*ones(betaA.n_elem,1) );
      jA = join_cols(jA, conv_to<mat>::from(A) ) ;
      if (intercept == 1) mu[m] = dot(betaA, xbar.elem(A)) ;
    }
  }
  if (!naive) {
    nonzeros *= 1+lambda2;
    mu = ybar - (1+lambda2) * mu;
  } else {
    mu = ybar - mu;
  }
  
  // Updating monitored quantities
  if (monitor > 0) D_star = J_hat - J_star;
  
  return List::create(
    Named("nzeros")     = nonzeros ,
    Named("iA")         = iA       ,
    Named("jA")         = jA       ,
    Named("mu")         = mu       ,
    Named("lambda1")    = lambda1  ,
    Named("df")         = df       ,
    Named("monitoring") = 
      List::create(
        Named("it.active")      = it_active,
        Named("it.optim")       = it_optim ,
        Named("max.grd")        = max_grd  ,
        Named("converge")       = converge ,
        Named("pensteps.timer") = timing, 
        Named("delta.hat")      = D_hat    ,
        Named("delta.star")     = D_star   
      )
  );
  
}
