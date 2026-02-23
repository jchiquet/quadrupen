/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _ElasticNet_H
#define _ElasticNet_H

#include "GenericRegularizer.h"
#include "OptimizerL1.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class ElasticNet : 
  public GenericRegularizer<matrix,Norm::L1>{
  public:
    
    using GenericRegularizer<matrix,Norm::L1>::intercept_ ;
    using GenericRegularizer<matrix,Norm::L1>::lambdas_   ;
    using GenericRegularizer<matrix,Norm::L1>::penalty_   ;
    using GenericRegularizer<matrix,Norm::L1>::set_  ;
    using GenericRegularizer<matrix,Norm::L1>::data_ ;
    using GenericRegularizer<matrix,Norm::L1>::df_   ;
    using GenericRegularizer<matrix,Norm::L1>::lambda_factor_ ;
    using GenericRegularizer<matrix,Norm::L1>::get_lambda_seq ;
  
  ElasticNet(const RegressionData<matrix>&, const List&, const List&);
  
  List solution_path(const List&);

  void optimality_gap(double lambda_, uword type) ;
  
  const sp_mat coefficients() const { 
    return sp_mat(join_cols(iA_, jA_), nzeros_, lambdas_.size(), data_.p_) ; 
  }

  const sp_mat debiased_coefficients() const { 
    return sp_mat(join_cols(iA_, jA_), debiased_, lambdas_.size(), data_.p_) ; 
  }

  const vector<double>& intercept_debiased() const { return intercept_debiased_ ; }
  
  const sp_mat active_var() const { 
    return sp_mat(join_cols(iA_, jA_), vec(iA_.n_elem, fill::ones), lambdas_.size(), data_.p_) ; 
  }
  
  // Specific to Elastic-Net regularization
  OptimizerL1<matrix> solver_ ; // Solvers for L1 penalty
  double gamma_   ; // overall amount of l2 penalty
  vec beta_       ; // vector of current parameters
  vec grad_       ; // vector of current gradient (smooth part)
  double J_       ; // current optimality gap
  double D_       ; // current move in the optimality gap
  vec   nzeros_   ; // contains non-zero value of beta
  vec   debiased_ ; // contains the debiased non-zero value of beta
  vector<double >intercept_debiased_ ; // contains the debiased vector of intercept
  urowvec iA_     ; // contains row indices of the non-zero values
  urowvec jA_     ; // contains column indices of the non-zero values

  // Compute degrees of freedom for the current estimate
  double get_df() ;

};

template <typename matrix>
ElasticNet<matrix>::ElasticNet(
  const RegressionData<matrix>& data, const List& regParam, const List& control) :
  GenericRegularizer<matrix,Norm::L1>::GenericRegularizer(data, regParam) {
    
    // set the penalty to l1
    penalty_ = Penalty<Norm::L1>() ;
    lambda_factor_ = as<vec>(regParam["lambda_factor"]) ;
    get_lambda_seq(regParam) ;

    // Set up the optimizer
    solver_ = OptimizerL1<matrix>(penalty_) ;
    
    // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
    gamma_   = as<double>(regParam["gamma"]) ;
    data_.scale_struct(sqrt(gamma_)*pow(lambda_factor_,-1)) ;
    data_.scale_regressors(lambda_factor_) ;
    
    // Initialize the active set, beta_ and gradient with starting coefficient
    vec beta0 = control["beta0"] ;
    uvec A0 = find(beta0) ;
    if (A0.is_empty()) {
      set_  = ActiveSet(data_, as<bool>(control["usechol"])) ;
      grad_ = - data_.XTy_ ;
    } else {
      set_  = ActiveSet(data_, A0, as<bool>(control["usechol"])) ;
      beta_ = beta0(A0) ;
      grad_ = - data_.XTy_ + set_.XTXA_ * beta_  ;
    }
  }

template <typename matrix>
double ElasticNet<matrix>::get_df() {
  
  double df = set_.size() + data_.centered_ ;
  if (gamma_ > 0) {
    // loop due to sparse encoding. should iterate over the n_zeros only...
    mat SAA(set_.size(),set_.size()) ;
    for (uword i=0;i<set_.size();i++){
      for (uword j=i;j<set_.size();j++){
        SAA(i,j) = data_.S_.at(set_.A_(i),set_.A_(j));
        SAA(j,i) = SAA(i,j);
      }
    }
    df -= trace(SAA * set_.XATXAinv_);
  }
  
  return(df);
}

template <typename matrix>
List ElasticNet<matrix>::solution_path(const List& control) {
  
  // Parameters controlling the optimization
  const bool verbose(control["verbose"])      ; // verbosity level
  const double accuracy(control["threshold"]) ; // precision required
  const uword maxiter(control["maxiter"])     ; // max # of passes in the active set
  const uword maxfeat(control["maxfeat"])     ; // max # of variables activated
  const uword monitoring(control["monitor"])  ; // optimality monitor (0=none; 1=Grandvalet; 2=Fenchel)
  
  SolverType algorithm = QUADRA; // Optimizer (default to QUADRA)
  if (as<std::string>(control["method"]) == "FISTA") algorithm = FISTA;

  // Variables monitoring the algorithm
  vector<double> gap, timing, J_hat, D_hat ; // timings and optimality measures
  vector<uword> status, iactive, ioptim    ; // convergence and # of inner/outer iterates

  // LAMBDA LOOP
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda_ : lambdas_) {
    if (verbose) {
      Rprintf("\n lambda_l1 = %f",lambda_) ;
      Rprintf("\n nb active variables = %i\n", set_.size()) ;
    }
    
    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    
    // dual norm of the gradient
    vec grd_norm = abs(grad_) - lambda_ ;
    grd_norm(set_.A_) = abs(grad_(set_.A_) + lambda_ * sign(beta_)) ;
    // variable associated with the highest violation of KKT conditions 
    uword var_in = grd_norm.index_max() ;
    double current_gap = std::max(0.0, grd_norm(var_in)) ;
    uvec zeroed ;
    
    uword current_it = 0 ; bool success = true ; 
    J_ = datum::inf ; D_ = datum::inf ;
    while ((current_gap > accuracy) && (current_it <= maxiter)) {
      R_CheckUserInterrupt();
      current_it++;

      // ________________________________________________________________________
      // VARIABLE ACTIVATION IF APPLICABLE
      //
      if (set_.is_in_[var_in] == 0) { // Is var_in already in the active set?
        set_.add_var(var_in, data_) ;
        beta_.resize(beta_.size()+1) ; // update the vector of active parameters
        beta_.tail(1) = - 1e-3 * sign(grad_(var_in)) ;
        if (verbose) {Rprintf("\tnewly added variable %i\n",var_in);}
      } else if (verbose) {Rprintf("\talready in %i\n",var_in);}
      
      // ________________________________________________________________________
      // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      //
      if (algorithm == FISTA) {
        ioptim.push_back(
          solver_.fista(beta_, lambda_, data_, set_, 1e-10, 10000)
        );
      } else { // QUADRA solver
        try {
          ioptim.push_back(
            solver_.quadratic_enet(beta_, lambda_, data_, set_, 1e-5, 10000)
          );
        } catch (std::runtime_error& error) {
          if (verbose > 0) {
            Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
          }
          success = false ;
        }
      }
      
      // OPTIMALITY TESTING
      grad_ = - data_.XTy_ + set_.XTXA_ * beta_ ;
      grd_norm = penalty_.elt_norm(grad_) - lambda_ ;
      grd_norm(set_.A_) = abs(grad_(set_.A_) + lambda_ * sign(beta_)) ;
      var_in = grd_norm.index_max() ;
      current_gap = std::max(0.0, grd_norm(var_in)) ;
      
      if (monitoring > 0) {
        optimality_gap(lambda_, monitoring) ;
        J_hat.push_back(J_) ;
        D_hat.push_back(D_) ;
      }
      
    }
    if (verbose) Rprintf("\tcurrent gap = %f\n",current_gap) ;

    // Checking convergence status
    gap.push_back(current_gap) ;
    iactive.push_back(current_it) ;
    status.push_back(0) ;
    if (current_it >= maxiter)       { status.back() = 1 ; }
    if (set_.size() > maxfeat) { status.back() = 2 ; }
    if (!success)                    { status.back() = 3 ; }
    
    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      set_.inverse_Gram() ;
      nzeros_   = join_cols(nzeros_, beta_/(data_.norm_X_(set_.A_) % lambda_factor_(set_.A_)));
      vec beta_debiased = set_.XATXAinv_ * (data_.XTy_(set_.A_) - data_.X_bar_(set_.A_) * accu(data_.y_)) ;
      debiased_ = join_cols(debiased_, beta_debiased/(data_.norm_X_(set_.A_) % lambda_factor_(set_.A_)));
      intercept_.push_back(data_.y_bar_ - dot(beta_, data_.X_bar_(set_.A_)));
      intercept_debiased_.push_back(data_.y_bar_ - dot(beta_debiased, data_.X_bar_(set_.A_))) ;
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
      Named("gap_hat")        = J_hat  ,
      Named("delta_hat")      = D_hat  ,
      Named("convergence")    = status ,
      Named("pensteps_timer") = timing
    )
  );
}


template <typename matrix>
void ElasticNet<matrix>::optimality_gap(double lambda, uword type) {

  // gamma equals the max |gradient|
  double nu = norm(grad_, "inf");
  double loss = .5 * pow(data_.norm_y_ ,2) + 
    dot(beta_, .5 * set_.XATXA_ * beta_ - data_.XTy_(set_.A_)) ;
  double old_J = J_, old_D = D_ ;
  J_ = loss - dot(beta_, grad_(set_.A_))  ;
  uvec Ac ;
  
  switch (type) {
  case 1: // Grandvalet's bound
    Ac = find(grad_ > nu); // set of adversarial variables outside the boundary
    D_ = J_ * (1 - lambda/nu) - 
      (pow(lambda,2)/(2*gamma_))*((lambda*(data_.p_-Ac.n_elem))/nu + 
      pow(norm(grad_(Ac),2)/nu,2)-data_.p_);
    break;
  case 2: // Fenchel's bound
    if (nu < lambda) nu = lambda;
    D_ = loss * (1+pow(lambda/nu,2)) + sum(abs(lambda*beta_)) + 
      (lambda/nu)*(dot(beta_,data_.XTy_(set_.A_))-pow(data_.norm_y_,2));
    break;
  default: 
    D_ = datum::inf ;
    break;
  }
  
  // keep the smallest bound reached so far for a given lambda value
  if ((old_J < J_) && (old_D - D_) < (old_J - J_)) {D_ = old_D ; }
  
}

#endif

