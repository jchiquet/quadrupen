/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "BoundedRegression.h"

using namespace Rcpp;
using namespace arma;

BoundedRegression::BoundedRegression(
  RegressionData<mat>& data, const List& regParam, const List& control) :
  Regularizer<mat>::Regularizer(data, regParam) {
    
    // set the penalty to l infinity
    penalty_ = SimplePenalty<SimpleNorm::LINF>() ;
    get_lambda_seq(get_lambda_max(), regParam) ;

    // Set up the optimizer 
    solver_ = OptimizerLINF<mat>(penalty_, control) ;
    
    // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
    data_.scale_struct(sqrt(gamma_)*ones(data_.p_)) ;
    // data_.scale_struct(sqrt(gamma_)*pow(lambda_factor_,-1/2)) ;
    
    // Initialize the active set with starting coefficient
    set_= ActiveSet(data, as<bool>(control["usechol"])) ;
    
    // Compute the Gram matrix (+ S scaled)
    data_.precompute_XTX() ;

    beta_ = zeros<vec>(data_.p_) ; // vector of current parameters
    grad_ = -data_.XTy_          ; // vector of current gradient (smooth part)
  }

double BoundedRegression::get_df() {

  uvec U = find(abs(beta_) < max(abs(beta_))) ;
  double df = data_.centered_ + U.size();
  
  mat SUU(U.size(), U.size()) ;
  if (gamma_ > 0) {
    mat C = inv_sympd(data_.XTX_(U,U));
    // loop due to sparse encoding. should iterate over the n_zeros only...
    for (uword i=0;i<U.size();i++){
      for (uword j=i;j<U.size();j++){
        SUU(i,j) = data_.S_.at(U(i),U(j));
        SUU(j,i) = SUU(i,j);
      }
    }
    df -= trace(SUU * C);
  }
  
  return(df);
}

List BoundedRegression::solution_path(const List& control) {

  // Parameters controlling the optimization
  const bool verbose(control["verbose"])      ; // verbosity level
  const double accuracy(control["threshold"]) ; // precision required
  const uword maxiter(control["maxiter"])     ; // max # of passes in the active set
  const uword maxfeat(control["maxfeat"])     ; // max # of variables activated

  SolverType algorithm = QUADRA; // Optimizer (default to QUADRA)
  if (as<std::string>(control["method"]) == "FISTA") {
    algorithm = FISTA;
    set_.reset();
    set_.add_vars(all, data_);
  }

  // Variables monitoring the algorithm
  vector<double> gap, timing ; // timings and optimality measures
  vector<uword> status, iactive, ioptim ; // convergence and # of inner/outer iterates

  auto prox = [this](vec x, double l) {
    return(penalty_.proximal(x, l, this->lambda_factor_));
  } ;

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
        ioptim.push_back(
          solver_.fista_LM(beta_, grad_, lambda_, data_, set_, prox, 1e-3, 10000)
        );
        break;
      } else { // QUADRA solver
        try {
          if (current_it == maxiter) {
            throw std::runtime_error("Fail to converge...");
          } else {
            ioptim.push_back(
              solver_.quadratic(beta_, grad_, lambda_, data_, set_, accuracy, 10000)
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
      grad_ = - data_.XTy_ + data_.XTX_ * beta_ ;
      current_gap = penalty_.dual_norm(grad_, lambda_factor_) - lambda_ ;
    } while ((current_gap > accuracy) && (current_it <= maxiter));

    // Checking convergence status
    gap.push_back(fmax(0.0, penalty_.dual_norm(grad_, lambda_factor_) - lambda_)) ;
    iactive.push_back(current_it) ;
    status.push_back(0) ;
    if (current_it >= maxiter) { status.back() = 1 ; }
    if ((set_.size() > maxfeat) & 
        (algorithm == QUADRA)) { status.back() = 2 ; }

    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      coef_ = join_rows(coef_, beta_/data_.norm_X_) ;
      intercept_.push_back(data_.y_bar_ - as_scalar(dot(beta_, data_.X_bar_)));
      df_.push_back(get_df()) ;
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
