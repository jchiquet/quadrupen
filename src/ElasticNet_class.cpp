/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "ElasticNet_class.h"

using namespace Rcpp;
using namespace arma;

ElasticNet::ElasticNet(
  RegressionData<mat>& data, const bool& intercept, const List& regParam, const List& control) :
  GenericRegularizer<mat>::GenericRegularizer(data, intercept, regParam) {
    
    // set the penalty to l1
    penalty_ = Penalty(L1) ;
    lambda_factor_ = as<vec>(regParam["lambda_factor"]) ;
    get_lambda_seq(regParam) ;
    
    // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
    gamma_   = as<double>(regParam["gamma"]) ;
    data_.scale_struct(sqrt(gamma_)*pow(lambda_factor_,-1/2)) ;

    // Initialize the active set and beta_ with starting coefficient
    // vec beta0 = control["beta0"] ;
    // set_= ActiveSet(data, find(abs(beta0)), as<bool>(control["usechol"])) ;
    // beta_ = beta0.elem(set_.A_) ;

    set_= ActiveSet(data, as<bool>(control["usechol"])) ;

    // smooth part of the gradient
    grad_ = -data_.XTy_ ;
    if (!set_.A_.is_empty()) grad_ += set_.XTXA_ * beta_  ;

  }

double ElasticNet::get_df() {
  
  mat SAA(set_.size(),set_.size()) ;
  double df = set_.size() ;
  
  if (gamma_ > 0) {
    mat B ;
    if (!set_.R_.is_empty()) {
      B = solve(trimatu(set_.R_), eye(set_.R_.n_cols, set_.R_.n_cols));
      B = B * B.t();
    } else {
      B = inv_sympd(set_.XATXA_);
    }
    // loop due to sparse encoding. should iterate over the n_zeros only...
    for (uword i=0;i<set_.size();i++){
      for (uword j=i;j<set_.size();j++){
        SAA(i,j) = data_.S_.at(set_.A_(i),set_.A_(j));
        SAA(j,i) = SAA(i,j);
      }
    }
    df -= trace(SAA * B);
  }
  
  return(df);
}

List ElasticNet::solution_path(const List& control) {

  // Parameters controlling the optimization
  const bool verbose(as<bool>(control["verbose"]))        ; // verbosity level
  const double accuracy(as<double>(control["threshold"])) ; // precision required
  const uword maxiter(as<uword>(control["maxiter"]))      ; // max # of passes in the active set
  const uword maxfeat(as<uword>(control["maxfeat"]))      ; // max # of variables activated

  SolverType algorithm = QUADRA; // Optimizer (default to QUADRA)
  if (as<std::string>(control["method"]) == "FISTA") {
    algorithm = FISTA;
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
     if (verbose) {
       Rprintf("\n lambda_l1 = %f",lambda_) ;
       Rprintf("\n nb active variables = %i\n",set_.active().n_elem) ;
     }

    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE

    // dual norm of gradient for inactive variables
    vec grd_norm = penalty_.elt_norm(grad_) - lambda_ ;
    // gradient for active variables
    grd_norm.elem(set_.A_) = abs(grad_(set_.A_) + lambda_ * sign(beta_)) ;
    double current_gap = std::max(0.0, grd_norm.max()) ;

    uvec  null ; // stores the variables which go to zero during optimization
    uword current_it = 0 ;
    bool success = true ;
    while ((current_gap > accuracy) && (current_it <= maxiter)) {
      R_CheckUserInterrupt();
      current_it++;

      // ________________________________________________________________________
      // VARIABLE ACTIVATION IF APPLICABLE
      //
      // variable associated with the highest violation of KKT conditions 
      uword var_in = grd_norm.index_max() ;
      if (set_.is_in_[var_in] == 0) { // Is var_in already in the active set?
        set_.add_var(var_in, data_) ;
        beta_.resize(beta_.size()+1) ; // update the vector of active parameters
        beta_.tail(1) = 0.0   ;
        if (verbose) {Rprintf("\tnewly added variable %i\n",var_in);}
      }

      // ________________________________________________________________________
      // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      //
      if (algorithm == FISTA) {
        Optimizer solver(data_, penalty_, algorithm) ;
        ioptim.push_back(
          solver.run(beta_, lambda_, set_, set_.XATXA_, 1e-10, 10000)
        );
      } else { // QUADRA solver
        try {
          Optimizer solver(data_, penalty_, algorithm) ;
          ioptim.push_back(
            solver.run(beta_, lambda_, set_, set_.XATXA_, 1e-10, 50)
          );
        } catch (std::runtime_error& error) {
          if (verbose > 0) {
            Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
          }
          success = false ;
        }
      }
      
      // VARIABLE DELETION IF APPLICABLE
      
      // update the smooth part of the gradient
      grad_ = - data_.XTy_ + set_.XTXA_ * beta_ ;
      // removing variables zeroed during optimization
      null = sort(find(abs(beta_) + abs(grad_(set_.A_)) - lambda_ < ZERO), "descend") ;
      if (!null.is_empty()) {
        set_.del_vars(null) ;
        for (uword j=0; j<null.n_elem; j++) {
          if (verbose) Rprintf("\tremoving variable %i\n",null[j]);
          beta_.shed_row(null[j]) ;
        }
      }

      // OPTIMALITY TESTING
      // dual norm of gradient for inactive variables
      grd_norm = penalty_.elt_norm(grad_) - lambda_ ;
      // gradient for active variables
      grd_norm.elem(set_.A_) = abs(grad_(set_.A_) + lambda_ * sign(beta_)) ;
      current_gap = std::max(0.0, grd_norm.max()) ;
    }
    Rprintf("\tcurrent gap = %f\n",current_gap) ;
    
    // Checking convergence status
    gap.push_back(current_gap) ;
    iactive.push_back(current_it) ;
    status.push_back(0) ;
    if (current_it >= maxiter)  { status.back() = 1 ; }
    if (set_.size()  > maxfeat) { status.back() = 2 ; }
    if (!success)               { status.back() = 3 ; }

    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      nzeros_ = join_cols(nzeros_, beta_/(data_.norm_X_.elem(set_.A_) % lambda_factor_.elem(set_.A_)));
      const_.push_back(as_scalar(dot(beta_, data_.X_bar_.elem(set_.A_))));
      iA_ = join_rows(iA_, df_.size()*ones<urowvec>(set_.size()) );
      jA_ = join_rows(jA_, set_.A_.t()) ;
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
