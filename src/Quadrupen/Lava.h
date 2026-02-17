/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _Lava_H
#define _Lava_H

#include "GenericRegularizer.h"
#include "OptimizerL1.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class Lava : 
  public GenericRegularizer<matrix,Norm::L1>{
  
  using GenericRegularizer<matrix,Norm::L1>::lambdas_ ;
  using GenericRegularizer<matrix,Norm::L1>::penalty_ ;
  using GenericRegularizer<matrix,Norm::L1>::set_     ;
  using GenericRegularizer<matrix,Norm::L1>::data_    ;
  using GenericRegularizer<matrix,Norm::L1>::const_   ;
  using GenericRegularizer<matrix,Norm::L1>::df_      ;
  using GenericRegularizer<matrix,Norm::L1>::lambda_factor_ ;
  using GenericRegularizer<matrix,Norm::L1>::get_lambda_seq ;
  
public:
  
  Lava(const RegressionData<matrix>&, const RegressionData<matrix>&, const bool&, const List&, const List&);

  double get_lambda_max() {
    return(penalty_.dual_norm(scaled_data_.XTy_));
  } ;
  
  List solution_path(const List&);

  const double& struct_tuning() const { return gamma_ ; }
  
  const sp_mat sparse_coefficients() const { 
    return sp_mat(join_cols(iA_, jA_), nzeros_, lambdas_.size(), data_.p_) ; 
  }

  const vector<vec>& dense_coefficients() const { 
    return dense_coef_ ; 
  }
  
  const sp_mat active_var() const { 
    return sp_mat(join_cols(iA_, jA_), vec(iA_.n_elem, fill::ones), lambdas_.size(), data_.p_) ; 
  }
  
  // Specific to Elastic-Net regularization
  OptimizerL1<matrix> solver_ ; // Solvers for L1 penalty
  double gamma_ ; // overall amount of l2 penalty
  RegressionData<matrix> scaled_data_   ; // data structure
  vector<vec> dense_coef_ ; // matrix of dense coefficients
  vec beta_     ; // vector of current dense parameters
  vec delta_    ; // vector of current sparse parameters
  vec grad_     ; // vector of current gradient (smooth part)
  mat Proj_     ; // Lava projector matrix
  mat XTXinv_   ; // Lava projector matrix
  vec   nzeros_ ; // contains non-zero value of delta
  urowvec iA_   ; // contains row indices of the non-zero values
  urowvec jA_   ; // contains column indices of the non-zero values
  
  // Compute degrees of freedom for the current estimate
  double get_df() ;
  
};

template <typename matrix>
Lava<matrix>::Lava(
  const RegressionData<matrix>& data, 
  const RegressionData<matrix>& scaled_data, const bool& intercept, const List& regParam, const List& control) :
  GenericRegularizer<matrix,Norm::L1>::GenericRegularizer(data, intercept, regParam), 
  scaled_data_(scaled_data) {
    
    // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
    lambda_factor_ = as<vec>(regParam["lambda_factor"]) ;
    gamma_ = as<double>(regParam["gamma"]) ;

    // Compute the Projection matrices 
    mat XTX = data_.X_.t() * data_.X_ - data_.n_ * data_.X_bar_ * data_.X_bar_.t() + data_.S_;
    XTXinv_ = inv_sympd(XTX) ;
    Proj_ = data_.X_ * XTXinv_ * data_.X_.t() ;

    // set the penalty to L1 and initialize the optimizer
    penalty_ = Penalty<Norm::L1>() ;
    get_lambda_seq(regParam) ;
    solver_ = OptimizerL1<matrix>(penalty_) ; 
    
    // Initialize the active set, beta_ and gradient with starting coefficient
    vec beta0 = control["beta0"] ;
    uvec A0 = find(beta0) ;
    if (A0.is_empty()) {
      set_  = ActiveSet(scaled_data_, as<bool>(control["usechol"])) ;
      grad_ = - scaled_data_.XTy_ ;
    } else {
      set_  = ActiveSet(scaled_data_, A0, as<bool>(control["usechol"])) ;
      beta_ = beta0(A0) ;
      grad_ = - scaled_data_.XTy_ + set_.XTXA_ * beta_  ;
    }
  }

template <typename matrix>
double Lava<matrix>::get_df() {
  
  double df = set_.size() ;
  mat B ;
  if (set_.use_chol_) {
    B = set_.Rinv_ * set_.Rinv_.t();
  } else {
    B = inv_sympd(set_.XATXA_);
  }
  mat K = diagmat(ones(scaled_data_.n_)) - 
    scaled_data_.X_.cols(set_.A_) * B * scaled_data_.X_.cols(set_.A_).t() ;
  
  df -= trace(K * Proj_);
  
  return(df);
}

template <typename matrix>
List Lava<matrix>::solution_path(const List& control) {
  
  // Parameters controlling the optimization
  const bool verbose(control["verbose"])      ; // verbosity level
  const double accuracy(control["threshold"]) ; // precision required
  const uword maxiter(control["maxiter"])     ; // max # of passes in the active set
  const uword maxfeat(control["maxfeat"])     ; // max # of variables activated
  const uword monitoring(control["monitor"])  ; // optimality monitor (0=none; 1=Grandvalet; 2=Fenchel)
  
  SolverType algorithm = QUADRA; // Optimizer (default to QUADRA)
  if (as<std::string>(control["method"]) == "FISTA") algorithm = FISTA;
  
  // Variables monitoring the algorithm
  vector<double> gap, timing ; // timings and optimality measures
  vector<uword> status, iactive, ioptim    ; // convergence and # of inner/outer iterates
  
  // LAMBDA LOOP
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda_ : lambdas_) {
    if (verbose) {
      Rprintf("\n lambda_l1 = %f",lambda_) ;
      Rprintf("\n nb active variables = %i\n",set_.size()) ;
    }
    
    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    
    // dual norm of the gradient
    vec grd_norm = abs(grad_) - lambda_ ;
    grd_norm(set_.A_) = abs(grad_(set_.A_) + lambda_ * sign(delta_)) ;
    // variable associated with the highest violation of KKT conditions 
    uword var_in = grd_norm.index_max() ;
    double current_gap = std::max(0.0, grd_norm(var_in)) ;
    uvec zeroed ;
    
    uword current_it = 0 ; bool success = true ; 
    while ((current_gap > accuracy) && (current_it <= maxiter)) {
      R_CheckUserInterrupt();
      current_it++;
      
      // ________________________________________________________________________
      // VARIABLE ACTIVATION IF APPLICABLE
      //
      if (set_.is_in_[var_in] == 0) { // Is var_in already in the active set?
        set_.add_var(var_in, scaled_data_) ;
        delta_.resize(delta_.size()+1) ; // update the vector of active parameters
        delta_.tail(1) = - 1e-3 * sign(grad_(var_in)) ;
        if (verbose) {Rprintf("\tnewly added variable %i\n",var_in);}
      } else if (verbose) {Rprintf("\talready in %i\n",var_in);}
      
      // ________________________________________________________________________
      // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      //
      if (algorithm == FISTA) {
        ioptim.push_back(
          solver_.fista(delta_, lambda_, scaled_data_, set_, 1e-10, 10000)
        );
      } else { // QUADRA solver
        try {
          ioptim.push_back(
            solver_.quadratic_enet(delta_, lambda_, scaled_data_, set_, 1e-5, 10000)
          );
        } catch (std::runtime_error& error) {
          if (verbose > 0) {
            Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
          }
          success = false ;
        }
      }
      
      // OPTIMALITY TESTING
      grad_ = - scaled_data_.XTy_ + set_.XTXA_ * delta_ ;
      grd_norm = penalty_.elt_norm(grad_) - lambda_ ;
      grd_norm(set_.A_) = abs(grad_(set_.A_) + lambda_ * sign(delta_)) ;
      var_in = grd_norm.index_max() ;
      current_gap = std::max(0.0, grd_norm(var_in)) ;

    }
    if (verbose) Rprintf("\tcurrent gap = %f\n",current_gap) ;
    
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
      nzeros_ = join_cols(nzeros_, delta_/(scaled_data_.norm_X_(set_.A_) % lambda_factor_(set_.A_)));
      if (set_.size() > 0) {
        beta_ = (XTXinv_ * ( data_.XTy_ - data_.X_.t() * data_.X_.cols(set_.A_) * delta_)) / data_.norm_X_ ;
      } else {
        beta_ = (XTXinv_ * data_.XTy_) / data_.norm_X_ ;
      }
      dense_coef_.push_back(beta_) ;
      const_.push_back(data_.y_bar_ - dot(delta_, data_.X_bar_(set_.A_)) - dot(beta_, data_.X_bar_));
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

#endif
