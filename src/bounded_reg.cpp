/*
 * Author: Julien CHIQUET
 *         julien.chiquet@inrae.fr
 */

#include "quadrupen_headers.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List bounded_reg_cpp(
    const Environment &dataModel   , // data structure
    const List        &tuningParam , // List of tuning parameters
    const List        &control       // config of the optimisation 
  ) {
  
  // Reading input variables
  const uword n             = dataModel["n"]  ; // sample size
  const uword p             = dataModel["d"]  ; // problem size
  const SEXP &X             = dataModel["X"]  ; // design matrix
  const arma::vec &y        = dataModel["y"]  ; // response vector
  const arma::sp_mat& S     = dataModel["S"]  ; // Structuring matrix
  const arma::vec &penscale = dataModel["wx"] ;  // penalty weights
  const arma::vec &xty      = dataModel["xty"]    ; // responses to predictors vector
  const arma::vec xbar      = dataModel["mean_X"] ; // mean of the predictors
  const arma::vec &normx    = dataModel["norm_X"] ; // norm of the predictors
  const double ybar         = dataModel["mean_y"] ; // mean of the predictors
  const arma::vec& weights  = dataModel["wy"]     ; // observation weights (not use at the moment)
  const bool sparse         = dataModel["sparse_encoding"] ; // boolean for sparse mode

  arma::vec lambda_linf      = tuningParam["linf"] ; // vector of LInf penalties
  const double lambda_l2     = tuningParam["l2"]   ; // scalar for the amount L2 penalty
  
  // double eps           = control["threshold"] ; // precision required
  const double eps           = control["threshold"] ; // precision required
  const double eps2 = pow(eps, 2) ;
  const arma::uword max_iter = control["max.iter"]  ; // max # of iterates of the active set
  const arma::uword max_feat = control["max.feat"]  ; // max # of variables activated
        arma::uword fun      = control["method"]    ; // solver (0=quadra, 1=pathwise, 2=fista)
  const arma::uword verbose  = control["verbose"]   ; // int for verbose mode (0/1/2)
  const bool bullet          = control["bulletproof"];// use Cholesky decomposition or not

  // STRUCTURATING MATRIX (embed lambda_l2)
  sp_mat diag_S = spdiags(sqrt(lambda_l2)*pow(penscale,-1/2), ivec({0}), p, p) ;
  const sp_mat S_lambda_l2 = diag_S * S * diag_S  ; // sparsely encoded structuring matrix
  
  // Managing the data matrix in both cases of sparse or dense coding
  mat x        ;
  mat xt       ;
  mat xtx      ;
  sp_mat sp_x  ;
  sp_mat sp_xt ;
  if (sparse) { // Check how x is encoded for reading
    sp_x = as<sp_mat>(X) ;
    sp_xt = sp_x.t() ;
    xtx = sp_xt * sp_x - n * xbar * xbar.t() ;
  } else {
    x = as<mat>(X) ;
    xt = x.t();
    xtx = xt * x - n * xbar * xbar.t() ;
  }
  xtx += S_lambda_l2 ;
  
  // Initializing "first level" variables (outside of the lambda_linf loop)
  uword n_lambda = lambda_linf.n_elem        ; // # of penalty levels
  colvec beta     = zeros<vec>(p)          ; // vector of current parameters
  uvec all(p);
  for (uword i=0;i<p;i++){all(i) = i;}
  uvec   B           = all                 ; // guys reaching the boundary
  urowvec iB                               ; // contains row indices of the bounded variables
  urowvec jB                               ; // contains column indices of the bounded variables
  uvec   I = setdiff(all,B)                ; // guys living in between the supremum
  mat    coef                              ; // matrice of solution path
  vec    grd                               ; // smooth part of the gradient
  vec    mu        = zeros<vec>(n_lambda)  ; // the intercept term
  vec    max_grd   = zeros<vec>(n_lambda)  ; // a vector with the successively reach duality gap
  vec    converge  = zeros<vec>(n_lambda)  ; // a vector indicating if convergence occured (0/1/2)
  uvec   it_active = zeros<uvec>(n_lambda) ; // # of loop in the active set for each lambda_linf
  uvec   it_optim                          ; // # of loop in the optimization process for each loop of the active se
  bool  success_optim = true               ; // was the internal system resolution successful?
  vec    timing      (n_lambda)            ; // succesive timing in
  vec  df          (n_lambda)              ; // degrees of freedom
  wall_clock timer                         ; // clock

  // Initializing "second level" variables (within the active set - for a fixed value of lamdba)
  int  nbr_opt  = 0                        ; // # of current calls to the optimization routine
  double L      = max (xtx.diag())         ; // Lipschitz coefficients
  
  // _____________________________________________________________
  //
  // START THE LOOP OVER LAMBDA
  timer.tic();
  for (uword m=0; m<n_lambda; m++) {
    if (verbose == 2) {Rprintf("\n lambda_linf = %f",lambda_linf(m)) ;}
    
    // _____________________________________________________________
    //
    // START THE ALGORITHM
    // _____________________________________________________________
    //
    
    // smooth part of the gradient
    grd     = -xty + xtx * beta        ;
    // dual norm of the gradient
    max_grd[m] = as_scalar(sum(abs(grd)) - lambda_linf[m]);
    if (max_grd[m] < 0) {
      max_grd[m] = 0;
    }
    
    while ((max_grd[m] > eps) && (it_active[m] <= max_iter)) {
      // _____________________________________________________________
      //
      // (1) KKT/SYSTEM RESOLUTION
      // _____________________________________________________________
      
      // save the number of iterates performed along the optimization process
      it_optim.reshape(nbr_opt + 1,1) ;
      
      switch (fun) {
      case 1 :
        it_optim[nbr_opt] = fista_breg(beta, xtx, xty, grd, lambda_linf[m], L, eps2);
        // Evaluating the set of variable reaching the boundary
        B    = find(abs(abs(beta) - max(abs(beta))) < ZERO );
        break;
      default:
        try {
          // If no convergence up to now...
          if (it_active[m] == max_iter) {
            throw std::runtime_error("Fail to converge...");
          } else {
            it_optim[nbr_opt] = quadra_breg(beta, xtx, xty, lambda_linf[m], grd, B, I);
          }
        } catch (std::runtime_error& error) {
          if (verbose > 0) {
            Rprintf("\nWarning: experiencing numerical instability... ");
          }
          if (bullet) {
            if (verbose > 0) {
              Rprintf("\nEntering 'bulletproof' mode: switching to proximal algorithm (slower but safer).");
            }
            it_active[m] = 0; // start this lambda all the way back
            fun = 1 ;
            it_optim[nbr_opt] = fista_breg(beta, xtx, xty, grd, lambda_linf[m], L, 1e-2);
            // reformating the output
            B = find(abs(abs(beta) - max(abs(beta))) < ZERO );
            I = setdiff(all, B) ;
          } else {
            if (verbose > 0) {
              Rprintf("\nCutting the solution path to this point, as you specified bulletproof=FALSE.");
            }
            success_optim = false ;
          }
        }
      }
      nbr_opt++;
      
      // _____________________________________________________________
      //
      // (2) OPTIMALITY TESTING
      // _____________________________________________________________
      
      // dual norm of gradient for unactive variable
      max_grd[m] = as_scalar(sum(abs(grd)) - lambda_linf[m]) ;
      if (max_grd[m] < 0) {
        max_grd[m] = 0;
      }
      
      // Moving to the next iterate
      it_active[m]++;
      
      // Cutting the path here if fail to converge or
      if (!success_optim) {
        break;
      }
      
      R_CheckUserInterrupt();
    }
    
    // Record the time ellapsed
    timing[m] = timer.toc() ;
    
    // Checking convergence status
    if (it_active[m] >= max_iter) {
      converge[m] = 1;
    }
    if (p-B.n_elem > max_feat) {
      converge[m] = 2 ;
    }
    if (!success_optim) {
      converge[m] = 3;
    }
    
    // degress of freedom
    df[m] = get_df_breg(lambda_l2, xtx, S_lambda_l2, I);
    
    // Preparing next value of the penalty
    if (converge[m] == 2 || converge[m] == 3) {
      lambda_linf = lambda_linf.subvec(0,m-1) ;
      converge    = converge.subvec(0,m)      ;
      max_grd     = max_grd.subvec(0,m-1)     ;
      it_active   = it_active.subvec(0,m)     ;
      timing      = timing.subvec(0,m)        ;
      df          = df.subvec(0,m)            ;
      break;
    } else {
      coef = join_rows(coef, beta/(normx % penscale));
      mu[m] = dot(beta, xbar) ;
      iB = join_rows(iB, m*ones<urowvec>(B.n_elem) );
      jB = join_rows(jB, B.t()) ;
    }
    
  }
  // END OF THE LOOP OVER LAMBDA
  
  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = lambda_linf,
      Named("l2") = lambda_l2 
    ),
    Named("beta")        = strans(coef),
    Named("active")      = sp_mat(join_cols(iB, jB), vec(iB.n_elem, fill::ones), lambda_linf.n_elem, p),
    Named("mu")          = ybar - mu   ,
    Named("df")          = df          ,
    Named("monitoring")  = 
      List::create(
        Named("it.active")      = it_active,
        Named("it.optim")       = it_optim ,
        Named("max.grd")        = max_grd  ,
        Named("converge")       = converge ,
        Named("pensteps.timer") = timing 
    )
  );
}


// ============================
// OLD ONE
// 

// // [[Rcpp::export]]
// Rcpp::List bounded_reg_old_cpp(SEXP X        ,
//                                const arma::vec& Y        ,
//                                const arma::sp_mat Struct  , // Structuring matrix
//                                SEXP LAMBDA1  ,
//                                const arma::uword N_LAMBDA ,
//                                const double MIN_RATIO,
//                                const arma::vec& PENSCALE ,
//                                double LAMBDA2  ,
//                                bool INTERCEPT,
//                                bool NORMALIZE,
//                                const arma::vec& WEIGHTS  ,
//                                bool NAIVE    ,
//                                double EPS      ,
//                                const arma::uword& MAXITER  ,
//                                const arma::uword& MAXFEAT  ,
//                                const arma::uword& FUN      ,
//                                int VERBOSE  ,
//                                bool SPARSE   ,
//                                bool BULLETPROOF) {
//   
//   // Reading input variables
//   bool intercept(INTERCEPT)   ; // boolean for intercept mode
//   bool normalize(NORMALIZE)   ; // boolean for standardizing the predictor
//   double lambda2(LAMBDA2)     ; // penalty levels
//   vec    weights(WEIGHTS)     ; // observation weights (not use at the moment)
//   vec    penscale(PENSCALE)    ; // penalty weights
//   bool   naive(NAIVE)       ; // naive elastic-net or not
//   vec    y(Y)           ; // reponse vector
//   double eps(EPS)         ; // precision required
//   uword  fun(FUN)         ; // solver (0=quadra, 1=pathwise, 2=fista)
//   int    verbose(VERBOSE)     ; // int for verbose mode (0/1/2)
//   bool   sparse(SPARSE)      ; // boolean for sparse mode
//   bool   bullet(BULLETPROOF) ; // int for verbose mode (0/1/2)
//   uword  max_iter(MAXITER)     ; // max # of iterates of the active set
//   uword  max_feat(MAXFEAT)     ; // max # of variables activated
//   
//   vec    xty   ; // responses to predictors vector
//   vec    xbar  ; // mean of the predictors
//   vec    meanx ; // mean of the predictors (rescaled)
//   vec    normx ; // norm of the predictors
//   double normy ; // norm of the response
//   double ybar  ; // mean of the response
//   uword n      ; // sample size
//   uword p      ; // problem size
//   
//   // Managing the data matrix in both cases of sparse or dense coding
//   mat x        ;
//   mat xt       ;
//   sp_mat sp_x  ;
//   sp_mat sp_xt ;
//   mat xtx      ;
//   if (sparse == 1) { // Check how x is encoded for reading
//     sp_x = as<sp_mat>(X) ;
//     standardize(sp_x, y, intercept, normalize, penscale, xty, normx, normy, xbar, ybar) ;
//     sp_xt = sp_x.t() ;
//     n = sp_x.n_rows ;
//     p = sp_x.n_cols ;
//     xtx = sp_xt * sp_x - n * xbar * xbar.t()  ;
//   } else {
//     x = as<mat>(X) ;
//     standardize(x, y, intercept, normalize, penscale, xty, normx, normy, xbar, ybar) ;
//     x  = x - sqrt(weights) * trans(xbar) ;
//     xt = x.t();
//     n = x.n_rows ;
//     p = x.n_cols ;
//     xtx = xt * x ;
//   }
//   meanx = xbar % penscale % normx;
//   
//   // STRUCTURATING MATRIX
//   sp_mat diag_S = spdiags(sqrt(lambda2)*pow(penscale,-1/2), ivec({0}), p, p) ;
//   sp_mat S = diag_S * Struct * diag_S  ; // sparsely encoded structuring matrix
//   xtx += S ; // S is scaled by lambda2
//   
//   // VECTOR OF TUNING PARAMETER FOR THE L1-PENALTY
//   vec lambda1 = get_lambda1(LAMBDA1, N_LAMBDA, MIN_RATIO, sum(abs(xty)));
//   uword n_lambda = lambda1.n_elem  ; // # of penalty levels
//   
//   // Initializing "first level" variables (outside of the lambda1 loop)
//   colvec beta     = zeros<vec>(p)          ; // vector of current parameters
//   uvec all(p);
//   for (uword i=0;i<p;i++){all(i) = i;}
//   uvec   B           = all                 ; // guys reaching the boundary
//   uvec   I = setdiff(all,B)                ; // guys living in between the supremum
//   mat    coef                              ; // matrice of solution path
//   vec    grd                               ; // smooth part of the gradient
//   vec    mu        = zeros<vec>(n_lambda)  ; // the intercept term
//   vec    max_grd   = zeros<vec>(n_lambda)  ; // a vector with the successively reach duality gap
//   vec    converge  = zeros<vec>(n_lambda)  ; // a vector indicating if convergence occured (0/1/2)
//   uvec   it_active = zeros<uvec>(n_lambda) ; // # of loop in the active set for each lambda1
//   uvec   it_optim                          ; // # of loop in the optimization process for each loop of the active se
//   bool  success_optim = true               ; // was the internal system resolution successful?
//   vec    timing      (n_lambda)            ; // succesive timing in
//   vec  df          (n_lambda)              ; // degrees of freedom
//   wall_clock timer                         ; // clock
//   mat   iB                                 ; // contains row indices of the bounded variables
//   mat   jB                                 ; // contains column indices of the bounded variables
//   
//   // Initializing "second level" variables (within the active set - for a fixed value of lamdba)
//   int  nbr_opt  = 0                        ; // # of current calls to the optimization routine
//   double L      = max (xtx.diag())         ; // Lipschitz coefficients
//   
//   // _____________________________________________________________
//   //
//   // START THE LOOP OVER LAMBDA
//   timer.tic();
//   for (uword m=0; m<n_lambda; m++) {
//     if (verbose == 2) {Rprintf("\n lambda1 = %f",lambda1(m)) ;}
//     
//     // _____________________________________________________________
//     //
//     // START THE ALGORITHM
//     // _____________________________________________________________
//     //
//     
//     // smooth part of the gradient
//     grd     = -xty + xtx * beta        ;
//     // dual norm of the gradient
//     max_grd[m] = as_scalar(sum(abs(grd)) - lambda1[m]);
//     if (max_grd[m] < 0) {
//       max_grd[m] = 0;
//     }
//     
//     while ((max_grd[m] > eps) && (it_active[m] <= max_iter)) {
//       // _____________________________________________________________
//       //
//       // (1) KKT/SYSTEM RESOLUTION
//       // _____________________________________________________________
//       
//       // save the number of iterates performed along the optimization process
//       it_optim.reshape(nbr_opt + 1,1) ;
//       
//       switch (fun) {
//       case 1 :
//         it_optim[nbr_opt] = fista_breg(beta, xtx, xty, grd, lambda1[m], L, pow(eps,2));
//         // Evaluating the set of variable reaching the boundary
//         B    = find(abs(abs(beta) - max(abs(beta))) < ZERO );
//         break;
//       default:
//         try {
//           // If no convergence up to now...
//           if (it_active[m] == max_iter) {
//             throw std::runtime_error("Fail to converge...");
//           } else {
//             it_optim[nbr_opt] = quadra_breg(beta, xtx, xty, lambda1[m], grd, B, I);
//           }
//         } catch (std::runtime_error& error) {
//           if (verbose > 0) {
//             Rprintf("\nWarning: experiencing numerical instability... ");
//           }
//           if (bullet) {
//             if (verbose > 0) {
//               Rprintf("\nEntering 'bulletproof' mode: switching to proximal algorithm (slower but safer).");
//             }
//             it_active[m] = 0; // start this lambda all the way back
//             eps = 1e-2 ; // with the proximal settings
//             fun = 1 ;
//             it_optim[nbr_opt] = fista_breg(beta, xtx, xty, grd, lambda1[m], L, pow(eps,2));
//             // reformating the output
//             B = find(abs(abs(beta) - max(abs(beta))) < ZERO );
//             I = setdiff(all, B) ;
//           } else {
//             if (verbose > 0) {
//               Rprintf("\nCutting the solution path to this point, as you specified bulletproof=FALSE.");
//             }
//             success_optim = false ;
//           }
//         }
//       }
//       nbr_opt++;
//       
//       // _____________________________________________________________
//       //
//       // (2) OPTIMALITY TESTING
//       // _____________________________________________________________
//       
//       // dual norm of gradient for unactive variable
//       max_grd[m] = as_scalar(sum(abs(grd)) - lambda1[m]) ;
//       if (max_grd[m] < 0) {
//         max_grd[m] = 0;
//       }
//       
//       // Moving to the next iterate
//       it_active[m]++;
//       
//       // Cutting the path here if fail to converge or
//       if (!success_optim) {
//         break;
//       }
//       
//       R_CheckUserInterrupt();
//     }
//     
//     // Record the time ellapsed
//     timing[m] = timer.toc() ;
//     
//     // Checking convergence status
//     if (it_active[m] >= max_iter) {
//       converge[m] = 1;
//     }
//     if (p-B.n_elem > max_feat) {
//       converge[m] = 2 ;
//     }
//     if (!success_optim) {
//       converge[m] = 3;
//     }
//     
//     // degress of freedom
//     df[m] = get_df_breg(lambda2, xtx, S, I);
//     
//     // Preparing next value of the penalty
//     if (converge[m] == 2 || converge[m] == 3) {
//       lambda1     =    lambda1.subvec(0,m-1) ;
//       converge    =  converge.subvec(0,m)    ;
//       max_grd     =   max_grd.subvec(0,m-1)  ;
//       it_active   = it_active.subvec(0,m)    ;
//       timing      =    timing.subvec(0,m)    ;
//       df          =    df.subvec(0,m)        ;
//       break;
//     } else {
//       if (any(penscale != 1)) {
//         coef = join_rows(coef,beta/(normx % penscale));
//       } else {
//         coef = join_rows(coef,beta/normx);
//       }
//       iB = join_cols(iB, m*ones(B.n_elem,1) );
//       jB = join_cols(jB, conv_to<mat>::from(B) );
//       if (intercept == 1) {
//         mu[m] = dot(beta, xbar) ;
//       }
//     }
//     
//   }
//   // END OF THE LOOP OVER LAMBDA
//   
//   if (!naive) {
//     coef *= 1+lambda2;
//     mu = ybar - (1+lambda2) * mu;
//   } else {
//     mu = ybar - mu;
//   }
//   
//   return List::create(Named("coefficients") = strans(coef),
//                       Named("iB")           = iB          ,
//                       Named("jB")           = jB          ,
//                       Named("mu")           = mu          ,
//                       Named("normx")        = normx       ,
//                       Named("lambda1")      = lambda1     ,
//                       Named("df")           = df          ,
//                       Named("it.active")    = it_active   ,
//                       Named("it.optim")     = it_optim    ,
//                       Named("max.grd")      = max_grd     ,
//                       Named("timing")       = timing      ,
//                       Named("converge")     = converge    );
// }
