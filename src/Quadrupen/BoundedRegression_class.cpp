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
    get_lambda_seq(regParam) ;
    
    // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
    gamma_   = as<double>(regParam["gamma"]) ;
    data_.scale_struct(sqrt(gamma_)*pow(lambda_factor_,-1/2)) ;
    
    // Initialize the active set with starting coefficient
    set_= ActiveSet(data, as<bool>(control["usechol"])) ;
    
    // Compute the Gram matrix (+ S scaled)  
    XTX = data_.X_.t() * data_.X_ - data_.n_ * data_.X_bar_ * data_.X_bar_.t() + data_.S_;
    
    beta_ = zeros<vec>(data_.p_) ; // vector of current parameters
    grad_ = -data_.XTy_          ; // vector of current gradient (smooth part)
    
  }

double BoundedRegression::get_df() {

  double bound = max(abs(beta_)) ;
  uvec A = find(abs(beta_) >= bound) ;
  
  mat SAA(A.size(),A.size()) ;
  double df = A.size() ;
  
  if (gamma_ > 0) {
    mat C = inv_sympd(XTX(A,A));
    // loop due to sparse encoding. should iterate over the n_zeros only...
    for (uword i=0;i<A.size();i++){
      for (uword j=i;j<A.size();j++){
        SAA(i,j) = data_.S_.at(A(i),A(j));
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
  if (as<std::string>(control["method"]) == "FISTA") {
    algorithm = FISTA;
    set_.reset();
    set_.add_vars(all, data_);
  }

  // Variables monitoring the algorithm
  vector<double> gap    ; // a vector with the successively reach duality gap
  vector<uword> status  ; // a vector indicating if convergence occurred (0/1/2)
  vector<uword> iactive ; // # of loop in the active set for each lambda
  vector<uword> ioptim  ; // # of loop in the optimization process for each loop of the active set
  vector<double> timing ; // successive timing for solving for each lambda value
   
  // LAMBDA LOOP
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda_ : lambdas_) {
    if (verbose) {Rprintf("\n lambda_linf = %f",lambda_) ;}

    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    uword current_it = 0 ;
    double current_gap = datum::inf ;
    do {
      R_CheckUserInterrupt();
      current_it++;
      if (algorithm == FISTA) {
        Optimizer solver(data_, penalty_, algorithm) ;
        ioptim.push_back(
          solver.run(beta_, lambda_, set_, XTX, 1e-3, 10000)
        );
        break;
      } else { // QUADRA solver
        try {
          if (current_it == maxiter) {
            throw std::runtime_error("Fail to converge...");
          } else {
            Optimizer solver(data_, penalty_, algorithm) ;
            ioptim.push_back(
              solver.run(beta_, lambda_, set_, XTX, accuracy, 10000)
            );
          }
        } catch (std::runtime_error& error) {
          if (verbose) {
            Rprintf("\nNumerical instability: switching to proximal algorithm (slower but safer).");
          }
          current_it = 0; // start this lambda all the way back, with FISTA algorithm
          algorithm = FISTA ;
          set_.reset();
          set_.add_vars(all, data_);
        }
      }

      // OPTIMALITY TESTING
      grad_ = - data_.XTy_ + XTX * beta_ ;
      current_gap = penalty_.dual_norm(grad_) - lambda_ ;
    } while ((current_gap > accuracy) && (current_it <= maxiter));

    // Checking convergence status
    gap.push_back(fmax(0.0, penalty_.dual_norm(grad_) - lambda_)) ;
    iactive.push_back(current_it) ;
    status.push_back(0) ;
    if (current_it >= maxiter) { status.back() = 1 ; }
    if ((set_.size() > maxfeat) & 
        (algorithm == QUADRA)) { status.back() = 2 ; }

    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      coef_.push_back(beta_/(data_.norm_X() % lambda_factor_)) ;
      const_.push_back(as_scalar(data_.y_bar_ - dot(beta_, data_.X_bar_)));
      df_.push_back(get_df()) ;
      iA_ = join_rows(iA_, df_.size()*ones<urowvec>(set_.size()) );
      jA_ = join_rows(jA_, set_.A_.t()) ;
    }

    timing.push_back(timer.toc()) ;
  } // END OF THE LOOP OVER LAMBDA

  lambdas_.resize(df_.size()) ;
  
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
