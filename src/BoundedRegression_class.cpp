/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "BoundedRegression_class.h"

using namespace Rcpp;
using namespace arma;

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
  
  // Initialize the set of variables living on the boundaries (all)
  B = regspace<uvec>(0, data_.p_-1) ;
  
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
  all = regspace<uvec>(0, data_.p_-1) ;
  uvec I = setdiff(all, B);
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
  const uword maxiter(as<uword>(control["maxiter"]))      ; // max # of passes in the active set
  const uword maxfeat(as<uword>(control["maxfeat"]))      ; // max # of variables activated

  SolverType algorithm = QUADRA; // Optimizer (default to QUADRA)
  if (as<std::string>(control["method"]) == "FISTA") {algorithm = FISTA;}

  // Variables monitoring the algorithm
  vector<double> gap    ; // a vector with the successively reach duality gap
  vector<uword> status  ; // a vector indicating if convergence occurred (0/1/2)
  vector<uword> iactive ; // # of loop in the active set for each lambda
  vector<uword> ioptim  ; // # of loop in the optimization process for each loop of the active set
  vector<double> timing ; // successive timing for solving for each lambda value
   
  // LAMBDA LOOP
  vec current_beta = zeros<vec>(data_.p_) ; // vector of current parameters
  wall_clock timer ; timer.tic(); // clock
  for(auto current_lambda : lambda_) {
     if (verbose) {Rprintf("\n lambda_linf = %f",current_lambda) ;}

    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    uword current_it = 0 ;
    bool  success = true ;
    vec current_grad ;
    double current_gap = datum::inf ;
    
    do {
      R_CheckUserInterrupt();
      current_it++;
      
      if (algorithm == FISTA) {
        Optimizer solver(data_, penalty_, algorithm) ;
        ioptim.push_back(
          solver.run(current_beta, B, current_grad, XTX, current_lambda, 1e-3, 10000)
        );
        break;
      } else { // QUADRA solver
        try {
          if (current_it == maxiter) {
            throw std::runtime_error("Fail to converge...");
          } else {
            Optimizer solver(data_, penalty_, algorithm) ;
            ioptim.push_back(
              solver.run(current_beta, B, current_grad, XTX, current_lambda, 1e-10, 50)
            );
          }
        } catch (std::runtime_error& error) {
          if (verbose) {
            Rprintf("\nNumerical instability: switching to proximal algorithm (slower but safer).");
          }
          current_it = 0; // start this lambda all the way back, with FISTA algorithm
          algorithm = FISTA ;
        }
      }

      // OPTIMALITY TESTING
      current_gap = std::max(0.0, penalty_.dual_norm(current_grad) - current_lambda) ;
    } while ((current_gap > accuracy) && (current_it <= maxiter) && (success));

    // Checking convergence status
    gap.push_back(std::max(0.0, penalty_.dual_norm(current_grad) - current_lambda)) ;
    iactive.push_back(current_it) ;
    status.push_back(0) ;
    if (current_it >= maxiter)         { status.back() = 1 ; }
    if (data_.p_ - B.n_elem > maxfeat) { status.back() = 2 ; }
    if (!success)                      { status.back() = 3 ; }

    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      beta_.push_back(current_beta/(data_.norm_X() % linf_weights)) ;
      mu_.push_back(as_scalar(data_.y_bar_ - dot(current_beta, data_.X_bar_)));
      df_.push_back(get_df()) ;
      iB_ = join_rows(iB_, df_.size()*ones<urowvec>(B.n_elem) );
      jB_ = join_rows(jB_, B.t()) ;
    }

    timing.push_back(timer.toc()) ;
  } // END OF THE LOOP OVER LAMBDA

  lambda_.resize(df_.size()) ;
  
  return(
    List::create(
      Named("it_active")      = iactive,
      Named("it_optim")       = ioptim ,
      Named("max_grd")        = gap    ,
      Named("convergence")    = status ,
      Named("pensteps_timer") = timing 
    )
  );
}
