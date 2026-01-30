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
    lambda_seq(regParam) ;
    
    // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
    gamma_   = as<double>(regParam["gamma"]) ;
    data_.scale_struct(sqrt(gamma_)*pow(lambda_factor_,-1/2)) ;

    // Initialize the active set with starting coefficient
    vec beta0 = control["beta0"] ;
    set_= ActiveSet(data, find(abs(beta0)), control["usechol"]) ;
    
    // WARM START
    // A_ = find(beta0 != 0) ;
    // betaA = beta0.elem(A) ;
    // if (A.n_elem > 0) {
    //   if (sparse) {
    //     xtxA = mat(sp_xt * sp_x.col(0)) - n * xbar * (xbar.elem(A)).t();
    //   } else {
    //     xtxA = mat(xt * x.cols(A)) - n * xbar * (xbar.elem(A)).t();
    //   }
    //   if (lambda_l2 > 0) {
    //     for (uword i=0; i<A.n_elem;i++) {
    //       xtxA.col(i) = xtxA.col(i) + S_lambda_l2.col(A(i));
    //       are_in(A(i)) = 1;
    //     }
    //   }
    //   grd += xtxA * betaA    ;
    //   nbr_in = A.n_elem      ;
    //   xAtxA = xtxA.rows(A)   ;
    //   if ((fun == 0) & (usechol)) {
    //     R = chol(xAtxA) ;
    //   }
    //   if (fun == 1) {
    //     xtxw(nbr_in) = dot(xAtxA.col(nbr_in),betaA);
    //   }
    // }
    
    
  }

double ElasticNet::get_df() {
  
  mat SAA(set_.size(),set_.size()) ;
  double df = set_.size() ;
  
  if (gamma_ > 0) {
    mat B ;
    if (!R_.is_empty()) {
      B = solve(trimatu(R_), eye(R_.n_cols, R_.n_cols));
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
     if (verbose) {
       Rprintf("\n lambda_l1 = %f",current_lambda) ;
       Rprintf("\n nb active variables = %i",set_.active().n_elem) ;
     }

    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE

    // smooth part of the gradient
    vec current_grad = -data_.XTy_ + set_.XTXA_ * current_beta.elem(set_.A_) ;
    // dual norm of gradient for inactive variables
    vec grd_norm = penalty_.elt_norm(current_grad) - current_lambda ;
    // gradient for active variables
    grd_norm.elem(set_.A_) = abs(current_grad(set_.A_) + current_lambda * sign(current_beta.elem(set_.A_))) ;
    double current_gap = std::max(0.0, grd_norm.max()) ;
    
    uvec  null ; // stores the variables which go to zero during optimization
    uword current_it = 0 ;
    bool success ;
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
        if (verbose) {Rprintf("newly added variable %i\n",var_in);}
      } else if (verbose) Rprintf("already in %i\n",var_in);

      // ________________________________________________________________________
      // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      //
      if (algorithm == FISTA) {
        Optimizer solver(data_, penalty_, algorithm) ;
        ioptim.push_back(
          solver.fista(current_beta, current_lambda, set_, set_.XATXA_, 1e-3, 10000)
        );
        break;
      } else { // QUADRA solver
        try {
          Optimizer solver(data_, penalty_, algorithm) ;
          ioptim.push_back(
            solver.fista(current_beta, current_lambda, set_, set_.XATXA_, 1e-10, 50)
          );
        } catch (std::runtime_error& error) {
          if (verbose > 0) {
            Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
          }
          success = false ;
        }
      }
      // update the smooth part of the gradient - WILL PROBABLY BE REMOVED
      current_grad = - data_.XTy_ + set_.XTXA_ * current_beta.elem(set_.A_) ;

      // 
      // VARIABLE DELETION IF APPLICABLE
      //
      // removing variables zeroed during optimization
      // null = sort(find(abs(beta) + abs(grad) - lambda < ZERO), "descend") ;
      if (!null.is_empty()) {
        for (uword j=0; j<null.n_elem; j++) {
          if (verbose) Rprintf("removing variable %i\n",null[j]);
          set_.del_var(null[j]) ;
        }
      }

      // OPTIMALITY TESTING
      // smooth part of the gradient
      current_grad = - data_.XTy_ + set_.XTXA_ * current_beta.elem(set_.A_) ;
      // dual norm of gradient for inactive variables
      grd_norm = penalty_.elt_norm(current_grad) - current_lambda ;
      // gradient for active variables
      grd_norm.elem(set_.A_) = abs(current_grad(set_.A_) + current_lambda * sign(current_beta.elem(set_.A_))) ;
      current_gap = std::max(0.0, grd_norm.max()) ;
    }

    // Checking convergence status
    gap.push_back(std::max(0.0, penalty_.dual_norm(current_grad) - current_lambda)) ;
    iactive.push_back(current_it) ;
    status.push_back(0) ;
    if (current_it >= maxiter)          { status.back() = 1 ; }
    if (set_.active().n_elem > maxfeat) { status.back() = 2 ; }
    if (!success)                       { status.back() = 3 ; }

    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      beta_.push_back(current_beta.elem(set_.A_)/ (data_.norm_X_.elem(set_.A_) % lambda_factor_.elem(set_.A_)));
      mu_.push_back(as_scalar(dot(current_beta.elem(set_.A_), data_.X_bar_.elem(set_.A_))));
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
