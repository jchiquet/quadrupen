/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "BoundedRegression_class.h"

using namespace Rcpp;
using namespace arma;

arma::uvec setdiff_bounded(arma::uvec& x, arma::uvec& y) {
  
  std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  
  return arma::conv_to<arma::uvec>::from(out);
}

BoundedRegression::BoundedRegression(
  RegressionData& data, const bool& intercept, const List& regParam) :
  data_ (data), intercept_ (intercept)
{

  // Initialize the penalty function
  penalty_ = Penalty() ;
  penalty_.setPenalty("linf") ;
  // TODO: include the weights in the penalty function
  linf_weights = as<vec>(regParam["linf_weights"]) ;
  // Generate the sequence of lambda (amount of linfty penalty)
  // In the future, should be a method of penalty to account for the weigths
  lambda_seq(regParam) ;
  
  // Initialize the set of variables: either on or within the l_infty boundaries
  all = regspace<uvec>(0, data_.p_-1) ;
  B = all ; // guys reaching the boundaries
  
  // Scale the structuring matrix according to the l_inf weights and the amount of l2 penalty 
  l2_penalty   = as<double>(regParam["l2"]) ;
  sp_mat diag_S = spdiags(sqrt(l2_penalty)*pow(linf_weights,-1/2), ivec({0}), data_.p_, data_.p_) ;
  S_scaled = diag_S * data_.S_ * diag_S  ; // sparsely encoded structuring matrix

  // Compute the Gram matrix (+ S_scaled)  
  XTX = data_.X_.t() * data_.X_ - data_.n_ * data_.X_bar_ * data_.X_bar_.t() + S_scaled;
}

void BoundedRegression::lambda_seq(const List& regParam) {
  if (regParam["linf"] != R_NilValue) {
    lambda_  = as<vector<double>>(regParam["linf"]) ;
  } else {
    double lambda_max = penalty_.dual_norm(data_.XTy_) ;
    lambda_ = conv_to<vector<double>>::from(logspace(log10(lambda_max),
                 log10(as<double>(regParam["min_ratio"])*lambda_max),
                 as<uword>(regParam["n_lambda"])
    ));
  }
}

double BoundedRegression::get_df() {
  double df ;
  mat C     ;
  mat SII(I.n_elem,I.n_elem) ;
  
  df = I.n_elem;
  if (l2_penalty > 0) {
    C = inv_sympd(XTX.submat(I,I));
    // loop due to sparse encoding.. should iterate over the n_zeros only...
    for (uword i=0;i<I.n_elem;i++){
      for (uword j=i;j<I.n_elem;j++){
        SII(i,j) = S_scaled.at(I(i),I(j));
        SII(j,i) = SII(i,j);
      }
    }
    df -= trace(SII * C);
  }
  
  return(df);
}

List BoundedRegression::solution_path(const List& control) {

  // Parameters controlling the optimization
  const bool verbose(as<bool>(control["verbose"]))        ; // verbosity level
  const double accuracy(as<double>(control["threshold"])) ; // precision required
  const uword maxiter(as<uword>(control["maxiter"]))      ; // max # of iterates of the active set
  const uword maxfeat(as<uword>(control["maxfeat"]))      ; // max # of variables activated
  const bool bullet(as<bool>(control["bulletproof"]))     ; // use Cholesky decomposition or not
  
  available_optimizer optimizer = QUADRA; // Optimizer (default to QUADRA)
  if (as<std::string>(control["method"]) == "FISTA") {optimizer = FISTA;}

  vector<double> max_grd  ; // a vector with the successively reach duality gap
  vector<uword> converge  ; // a vector indicating if convergence occurred (0/1/2)
  vector<uword> it_active ; // # of loop in the active set for each lambda
  vector<uword> it_optim  ; // # of loop in the optimization process for each loop of the active set
  vector<double> timing ; // successive timing for solving for each lambda value
   
  // LAMBDA LOOP
  vec current_beta = zeros<vec>(data_.p_) ; // vector of current parameters
  wall_clock timer ; // clock
  timer.tic();
  for(auto current_lambda : lambda_) {
     if (verbose) {Rprintf("\n lambda_linf = %f",current_lambda) ;}
  
    // smooth part of the gradient
    vec current_grad = -data_.XTy_ + XTX * current_beta ;
    // dual norm of the gradient
    // double current_max_grd = std::max(0.0, as_scalar(sum(abs(current_grad)) - current_lambda)) ;
    double current_max_grd = std::max(0.0, penalty_.dual_norm(current_grad) - current_lambda) ;
    
    // FIX-LAMBDA LOOP: IDENTIFY THE ACTIVE SET AND SOLVE
    uword current_it = 0 ;
    bool  success = true ;
    while ((current_max_grd > accuracy) && (current_it <= maxiter)) {
      // _____________________________________________________________
      //
      // (1) KKT/SYSTEM RESOLUTION
      // _____________________________________________________________
  
      switch (optimizer) {
      case FISTA :
        it_optim.push_back(
          solve_fista(current_beta, current_grad, current_lambda)
        );
        // Evaluating the set of variable reaching the boundary
        B = find(abs(penalty_.elt_norm(current_beta) - penalty_.pen_norm(current_beta)) < ZERO );
        break;
      default: // quadratic solver
        try {
          // If no convergence up to now...
          if (current_it == maxiter) {
            throw std::runtime_error("Fail to converge...");
          } else {
            it_optim.push_back(
              solve_quadratic(current_beta, current_grad, current_lambda)
            );
          }
        } catch (std::runtime_error& error) {
          if (verbose) {
            Rprintf("\nWarning: experiencing numerical instability... ");
          }
          if (bullet) {
            if (verbose) {
              Rprintf("\nEntering 'bulletproof' mode: switching to proximal algorithm (slower but safer).");
            }
            current_it = 0; // start this lambda all the way back
            optimizer = FISTA ;
            it_optim.push_back(
              solve_fista(current_beta, current_grad, current_lambda)
            );
            // reformatting the output
            B = find(abs(penalty_.elt_norm(current_beta) - penalty_.pen_norm(current_beta)) < ZERO );
          } else {
            if (verbose) {
              Rprintf("\nCutting the solution path to this point, as you specified bulletproof=FALSE.");
            }
            success = false ;
          }
        }
      }

      // _____________________________________________________________
      //
      // (2) OPTIMALITY TESTING
      // _____________________________________________________________

      // dual norm of gradient for inactive variable
      current_max_grd = std::max(0.0, penalty_.dual_norm(current_grad) - current_lambda) ;
      // current_max_grd = std::max(0.0, as_scalar(sum(abs(current_grad)) - current_lambda)) ;

      // Moving to the next iterate
      current_it++;

      // Cutting the path here if fail to converge or
      if (!success) { break; }

      R_CheckUserInterrupt();
    } // END OF THE SOLVER LOOP

    // Checking convergence status
    max_grd.push_back(current_max_grd) ;
    it_active.push_back(current_it) ;
    converge.push_back(0) ;
    if (current_it >= maxiter)         { converge.back() = 1 ; }
    if (data_.p_ - B.n_elem > maxfeat) { converge.back() = 2 ; }
    if (!success)                      { converge.back() = 3 ; }

    // Preparing next value of the penalty
    if (converge.back() >= 2) {
      break;
    } else {
      beta_.push_back(current_beta/(data_.norm_X() % linf_weights)) ;
      mu_.push_back(as_scalar(data_.y_bar_ - dot(current_beta, data_.X_bar_)));
      df_.push_back(get_df()) ;
      iB_ = join_rows(iB_, df_.size()*ones<urowvec>(B.n_elem) );
      jB_ = join_rows(jB_, B.t()) ;
    }

    // Record the time elapsed
    timing.push_back(timer.toc()) ;
  } // END OF THE LOOP OVER LAMBDA

  lambda_.resize(df_.size()) ;
  
  return(
    List::create(
      Named("it_active")      = it_active,
      Named("it_optim")       = it_optim ,
      Named("max_grd")        = max_grd  ,
      Named("convergence")    = converge ,
      Named("pensteps_timer") = timing 
    )
  );
}

int BoundedRegression::solve_quadratic(
    vec& beta,
    vec& grad, 
    double& lambda, 
    const uword& max_iter) {
  
  static const double zero = 2e-16     ;
  uword p        = beta.n_elem    ; // size of the problem
  uword iter     = 0              ; // count the number of systems solved
  double bound ; //
  uvec all = regspace<uvec>(0,p-1) ;
  uvec toB           ; // guys reaching the boundary after optimization
  uvec toI           ; // guys leaving the boundary after optimization
  vec  theta = -sign(grad.elem(B)) ;
  
  vec XX_B   ;
  mat XX     ;
  mat XX_II  ;
  mat R      ;
  vec b      ;
  vec tmp    ;
  
  toB = B;
  
  while ((toB.n_elem > 0) & (iter < max_iter)) {
    
    iter++;
    //
    // SOLVE THE QUADRATIC PROBLEM
    //
    
    // Constructing the system (KKT)
    XX_B = XTX.cols(B) * theta;
    
    if (I.is_empty()) {
      XX = sum(theta % XX_B.elem(B),0);
      b  = (dot(theta, data_.XTy_.elem(B)) - lambda);
      beta.elem(B) = theta * (b/XX) ;
      // keep a trace of the current boundary
    } else {
      if (I.n_elem > 1) {
        XX_II = XTX.submat(I,I) ;
      } else {
        XX_II = XTX(I,I) ;
      }
      XX   = join_rows(
        join_cols(sum(theta % XX_B.elem(B),0),XX_B.elem(I)),
        join_cols(strans(XX_B.elem(I)), XX_II));
      
      vec b = zeros<vec>(I.n_elem + 1) ;
      b[0] = dot(theta, data_.XTy_.elem(B)) - lambda;
      b.subvec(1,b.n_elem-1) = data_.XTy_.elem(I) ;
      
      // Solving with Cholesky factorization...
      R = chol(XX) ;
      tmp = solve(trimatu(R), solve(trimatl(strans(R)),b)) ;
      beta.elem(B) = theta * tmp[0] ;
      beta.elem(I) = tmp.subvec(1,tmp.n_elem-1) ;
    }
    // keep a trace of the current boundary
    bound = max(abs(beta.elem(B)));
    //
    // VARIABLES REACHING THE BOUNDARY
    //
    toB = find(abs(beta) > bound);
    beta.elem(toB) = bound * sign(beta.elem(toB));
    B = unique(join_cols(B,toB));
    I = setdiff_bounded(all,B);
    theta = sign(beta.elem(B)); // sign of the guys reaching the supremum
  }
  
  grad = -data_.XTy_ + XTX * beta ;
  
  //
  // VARIABLE LEAVING THE BOUDARY
  //
  toI = find(abs(theta + sign(grad.elem(B))) > zero);
  if (!toI.is_empty()) {
    toI = B.elem(toI);
  }
  if (!toI.is_empty()) {
    B = setdiff_bounded(B,toI);
    I = setdiff_bounded(all,B);
  }
  // If everyone is leaving the boundary, that's an issue...
  if (B.is_empty()) {
    throw std::runtime_error("Too much unstability");
  }
  
  return(iter) ;
}

int BoundedRegression::solve_fista(
    vec& beta0,
    vec& grad,
    double& lambda,
    const double& eps, 
    const uword& max_iter) {

  colvec betak = beta0  ; // output vector
  colvec betal = beta0     ;
  uword iter = 0          ; // current iterate
  double delta = 2*eps  ; // change in beta
  double L0 = max( XTX.diag()) ; // Lipchitz constant
  
  double t0 = 1.0, tk;
  bool found=false;
  double f0, fk ;
  double l_num, l_den ;
  
  while ((delta > eps*eps) && (iter < max_iter)) {
    
    f0 = as_scalar(.5 * strans(betal) * XTX * betal - strans(data_.XTy_) * betal) ;
    grad = -data_.XTy_ + XTX * betal ;
    
    // Line search over L
    while(!found) {
      betak = penalty_.proximal(betal - grad/L0, lambda /L0);
      
      fk = as_scalar(.5 * strans(betak) * XTX * betak - strans(data_.XTy_) * betak) ;
      l_num = as_scalar(2 * (fk - f0 - dot(grad, betak-betal) ));
      l_den = as_scalar(pow(norm(betak-betal,2),2));
      
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
    betal = betak + (t0-1)/tk * ( betak - beta0 );
    
    // preparing next iterate
    delta = sqrt(l_num);
    beta0 = betak;
    t0 = tk;
    found = false;
    iter++;
    
    R_CheckUserInterrupt();
  }
  
  return(iter) ;
  
}

