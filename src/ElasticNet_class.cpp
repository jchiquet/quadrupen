/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "ElasticNet_class.h"

using namespace Rcpp;
using namespace arma;

ElasticNet::ElasticNet(
  RegressionData<mat>& data, const bool& intercept, const List& regParam, SEXP BETA0) :
  GenericRegularizer<mat>::GenericRegularizer(data, intercept, regParam, BETA0) {
    
    // set the penalty to l infinity
    penalty_ = Penalty(L1) ;
    lambda_factor_ = as<vec>(regParam["lambda_factor"]) ;
    lambda_seq(regParam) ;
    
    // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
    gamma_   = as<double>(regParam["gamma"]) ;
    data_.scale_struct(sqrt(gamma_)*pow(lambda_factor_,-1/2)) ;
    
    // mat  R_       ; // Cholesky decomposition of XAtXA
    // vec  xtxw_    ; // t(x_A) * x_A * beta(A)
    
    // // Compute the Gram matrix (+ S scaled)  

    // Initialize the set of variables living on the boundaries (all)

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

// List ElasticNet::solution_path(const List& control) {
// 
//   // Parameters controlling the optimization
//   const bool verbose(as<bool>(control["verbose"]))        ; // verbosity level
//   const double accuracy(as<double>(control["threshold"])) ; // precision required
//   const uword maxiter(as<uword>(control["maxiter"]))      ; // max # of passes in the active set
//   const uword maxfeat(as<uword>(control["maxfeat"]))      ; // max # of variables activated
// 
//   SolverType algorithm = QUADRA; // Optimizer (default to QUADRA)
//   if (as<std::string>(control["method"]) == "FISTA") {algorithm = FISTA;}
// 
//   // Variables monitoring the algorithm
//   vector<double> gap    ; // a vector with the successively reach duality gap
//   vector<uword> status  ; // a vector indicating if convergence occurred (0/1/2)
//   vector<uword> iactive ; // # of loop in the active set for each lambda
//   vector<uword> ioptim  ; // # of loop in the optimization process for each loop of the active set
//   vector<double> timing ; // successive timing for solving for each lambda value
//    
//   // LAMBDA LOOP
//   vec current_beta = zeros<vec>(data_.p_) ; // vector of current parameters
//   wall_clock timer ; timer.tic(); // clock
//   for(auto current_lambda : lambda_) {
//      if (verbose) {Rprintf("\n lambda_linf = %f",current_lambda) ;}
// 
//     // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
//     uword current_it = 0 ;
//     vec current_grad ;
//     double current_gap = datum::inf ;
//     
//     do {
//       R_CheckUserInterrupt();
//       current_it++;
//       
//       if (algorithm == FISTA) {
//         Optimizer solver(data_, penalty_, algorithm) ;
//         ioptim.push_back(
//           solver.run(current_beta, A, current_grad, XTX, current_lambda, 1e-3, 10000)
//         );
//         break;
//       } else { // QUADRA solver
//         try {
//           if (current_it == maxiter) {
//             throw std::runtime_error("Fail to converge...");
//           } else {
//             Optimizer solver(data_, penalty_, algorithm) ;
//             ioptim.push_back(
//               solver.run(current_beta, A, current_grad, XTX, current_lambda, 1e-10, 50)
//             );
//           }
//         } catch (std::runtime_error& error) {
//           if (verbose) {
//             Rprintf("\nNumerical instability: switching to proximal algorithm (slower but safer).");
//           }
//           current_it = 0; // start this lambda all the way back, with FISTA algorithm
//           algorithm = FISTA ;
//         }
//       }
// 
//       // OPTIMALITY TESTING
//       current_gap = std::max(0.0, penalty_.dual_norm(current_grad) - current_lambda) ;
//     } while ((current_gap > accuracy) && (current_it <= maxiter));
// 
//     // Checking convergence status
//     gap.push_back(std::max(0.0, penalty_.dual_norm(current_grad) - current_lambda)) ;
//     iactive.push_back(current_it) ;
//     status.push_back(0) ;
//     if (current_it >= maxiter)         { status.back() = 1 ; }
//     // if (data_.p_ - B.n_elem > maxfeat) { status.back() = 2 ; }
// 
//     // Preparing next value of the penalty
//     if (status.back() >= 2) {
//       break;
//     } else {
//       beta_.push_back(current_beta/(data_.norm_X() % lambda_factor_)) ;
//       mu_.push_back(as_scalar(data_.y_bar_ - dot(current_beta, data_.X_bar_)));
//       df_.push_back(get_df()) ;
//       iA_ = join_rows(iA_, df_.size()*ones<urowvec>(A.n_elem) );
//       jA_ = join_rows(jA_, A.t()) ;
//     }
// 
//     timing.push_back(timer.toc()) ;
//   } // END OF THE LOOP OVER LAMBDA
// 
//   lambda_.resize(df_.size()) ;
//   
//   return(
//     List::create(
//       Named("it_active")      = iactive,
//       Named("it_optim")       = ioptim ,
//       Named("max_grd")        = gap    ,
//       Named("convergence")    = status ,
//       Named("pensteps_timer") = timing 
//     )
//   );
// }
