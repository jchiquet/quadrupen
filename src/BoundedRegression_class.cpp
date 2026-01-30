/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "BoundedRegression_class.h"

using namespace Rcpp;
using namespace arma;

BoundedRegression::BoundedRegression(
  RegressionData<mat>& data, const bool& intercept, const List& regParam, const List& control) :
  GenericRegularizer<mat>::GenericRegularizer(data, intercept, regParam) {

  // set the penalty to l infinity
  penalty_ = Penalty(LINF) ;
  lambda_factor_ = as<vec>(regParam["lambda_factor"]) ;
  lambda_seq(regParam) ;
  
  // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
  gamma_   = as<double>(regParam["gamma"]) ;
  data_.scale_struct(sqrt(gamma_)*pow(lambda_factor_,-1/2)) ;

  // Initialize the active set with starting coefficient
  set_= ActiveSet(data, as<bool>(control["usechol"])) ;
  
  // Compute the Gram matrix (+ S scaled)  
  XTX = data_.X_.t() * data_.X_ - data_.n_ * data_.X_bar_ * data_.X_bar_.t() + data_.S_;
  
}

double BoundedRegression::get_df() {

  mat SAA(set_.size(),set_.size()) ;
  double df = set_.size() ;
  
  if (gamma_ > 0) {
    mat C = inv_sympd(XTX(set_.A_,set_.A_));
    // loop due to sparse encoding. should iterate over the n_zeros only...
    for (uword i=0;i<set_.size();i++){
      for (uword j=i;j<set_.size();j++){
        SAA(i,j) = data_.S_.at(set_.A_(i),set_.A_(j));
        SAA(j,i) = SAA(i,j);
      }
    }
    df -= trace(SAA * C);
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
  vec current_grad = -data_.XTy_          ; // vector of current gradient (smooth part)
  wall_clock timer ; timer.tic(); // clock
  for(auto current_lambda : lambda_) {
    if (verbose) {Rprintf("\n lambda_linf = %f",current_lambda) ;}

    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    uword current_it = 0 ;
    double current_gap = datum::inf ;
    do {
      R_CheckUserInterrupt();
      current_it++;

      if (algorithm == FISTA) {
        Optimizer solver(data_, penalty_, algorithm) ;
        ioptim.push_back(
          solver.run(current_beta, current_lambda, set_, XTX, 1e-3, 10000)
        );
        break;
      } else { // QUADRA solver
        try {
          if (current_it == maxiter) {
            throw std::runtime_error("Fail to converge...");
          } else {
            Optimizer solver(data_, penalty_, algorithm) ;
            ioptim.push_back(
              solver.run(current_beta, current_lambda, set_, XTX, 1e-10, 10)
            );
          }
        } catch (std::runtime_error& error) {
          if (verbose) {
            Rprintf("\nNumerical instability: switching to proximal algorithm (slower but safer).");
          }
          current_it = 0; // start this lambda all the way back, with FISTA algorithm
          algorithm = FISTA ;
          // A.reset() ;
        }
      }

      // OPTIMALITY TESTING
      current_grad = - data_.XTy_ + XTX * current_beta ;
      current_gap = penalty_.dual_norm(current_grad) - current_lambda ;
    } while ((current_gap > accuracy) && (current_it <= maxiter));

    // Checking convergence status
    gap.push_back(fmax(0.0, penalty_.dual_norm(current_grad) - current_lambda)) ;
    iactive.push_back(current_it) ;
    status.push_back(0) ;
    if (current_it >= maxiter) { status.back() = 1 ; }
    if (set_.size() > maxfeat)    { status.back() = 2 ; }

    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      beta_.push_back(current_beta/(data_.norm_X() % lambda_factor_)) ;
      mu_.push_back(as_scalar(data_.y_bar_ - dot(current_beta, data_.X_bar_)));
      df_.push_back(get_df()) ;
      iA_ = join_rows(iA_, df_.size()*ones<urowvec>(set_.size()) );
      jA_ = join_rows(jA_, set_.A_.t()) ;
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
