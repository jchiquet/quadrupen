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

  using GenericRegularizer<matrix,Norm::L1>::lambdas_ ;
  using GenericRegularizer<matrix,Norm::L1>::penalty_ ;
  using GenericRegularizer<matrix,Norm::L1>::set_     ;
  using GenericRegularizer<matrix,Norm::L1>::data_    ;
  using GenericRegularizer<matrix,Norm::L1>::const_   ;
  using GenericRegularizer<matrix,Norm::L1>::df_      ;
  using GenericRegularizer<matrix,Norm::L1>::lambda_factor_ ;
  using GenericRegularizer<matrix,Norm::L1>::get_lambda_seq ;
  
  public:
  
  ElasticNet(const RegressionData<matrix>&, const bool&, const List&, const List&);
  
  List solution_path(const List&);

  double get_loss() {
    double loss = pow(data_.norm_y_ ,2) + dot(beta_, set_.XATXA_ * beta_) - 
      2*dot(beta_, data_.XTy_(set_.A_)) ;
    return (.5 * loss) ;
  };
  
  void optimality_gap(double lambda_, uword type) ;
  
  const double& struct_tuning() const { return gamma_ ; }

  const sp_mat coefficients() const { 
    return sp_mat(join_cols(iA_, jA_), nzeros_, lambdas_.size(), data_.p_) ; 
  }

  const sp_mat active_var() const { 
    return sp_mat(join_cols(iA_, jA_), vec(iA_.n_elem, fill::ones), lambdas_.size(), data_.p_) ; 
  }
  
  // Specific to Elastic-Net regularization
  OptimizerL1<matrix> solver_ ; // Solvers for L1 penalty
  double gamma_ ; // overall amount of l2 penalty
  vec beta_     ; // vector of current parameters
  vec grad_     ; // vector of current gradient (smooth part)
  double loss_  ; // current quadratic loss
  double J_     ; // current optimality gap
  double D_     ; // current move in the optimality gap
  vec   nzeros_ ; // contains non-zero value of beta
  urowvec iA_   ; // contains row indices of the non-zero values
  urowvec jA_   ; // contains column indices of the non-zero values

  // Compute degrees of freedom for the current estimate
  double get_df() ;

};

template <typename matrix>
ElasticNet<matrix>::ElasticNet(
  const RegressionData<matrix>& data, const bool& intercept, const List& regParam, const List& control) :
  GenericRegularizer<matrix,Norm::L1>::GenericRegularizer(data, intercept, regParam) {
    
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
  
  mat SAA(set_.size(),set_.size()) ;
  double df = set_.size() ;
  
  if (gamma_ > 0) {
    mat B ;
    if (set_.use_chol_) {
      // B = solve(trimatu(set_.R_), eye(set_.R_.n_cols, set_.R_.n_cols));
      // B = B * B.t();
      B = set_.Rinv_ * set_.Rinv_.t();
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

template <typename matrix>
List ElasticNet<matrix>::solution_path(const List& control) {
  
  // Parameters controlling the optimization
  const bool verbose(as<bool>(control["verbose"]))        ; // verbosity level
  const double accuracy(as<double>(control["threshold"])) ; // precision required
  const uword maxiter(as<uword>(control["maxiter"]))      ; // max # of passes in the active set
  const uword maxfeat(as<uword>(control["maxfeat"]))      ; // max # of variables activated
  const uword monitoring(as<uword>(control["monitor"]))   ; // optimality monitor (0=none; 1=Grandvalet; 2=Fenchel)
  
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
  vector<double> J_hat  ; // successive value of optimality gap
  vector<double> D_hat  ; // successive delta in optimality gap
  
  // LAMBDA LOOP
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda_ : lambdas_) {
    if (verbose) {
      Rprintf("\n lambda_l1 = %f",lambda_) ;
      Rprintf("\n nb active variables = %i\n",set_.size()) ;
    }
    
    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE
    
    // dual norm of gradient for inactive variables
    vec grd_norm = penalty_.elt_norm(grad_) - lambda_ ;
    // gradient for active variables
    grd_norm(set_.A_) = abs(grad_(set_.A_) + lambda_ * sign(beta_)) ;
    // variable associated with the highest violation of KKT conditions 
    uword var_in = grd_norm.index_max() ;
    double current_gap = std::max(0.0, grd_norm(var_in)) ;
    
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
        beta_.tail(1) = 0.0   ;
        if (verbose) {Rprintf("\tnewly added variable %i\n",var_in);}
      } else if (verbose) {Rprintf("\talready in %i\n",var_in);}
      
      // ________________________________________________________________________
      // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      //
      uvec  null ; // stores the variables which go to zero during optimization
      if (algorithm == FISTA) {
        ioptim.push_back(
          solver_.fista(beta_, lambda_, data_, set_, 1e-10, 10000)
        );
        null = find(abs(beta_) + abs(grad_(set_.A_)) - lambda_ < ZERO) ;
      } else { // QUADRA solver
        try {
          ioptim.push_back(
            solver_.quadratic_enet(beta_, lambda_, data_, set_, 1e-5, 10000)
          );
          null = find(abs(beta_) < ZERO) ;
        } catch (std::runtime_error& error) {
          if (verbose > 0) {
            Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
          }
          success = false ;
        }
      }
      
      // VARIABLE DELETION IF APPLICABLE
      if (!null.is_empty()) {
        null = sort(null, "descend") ;
        set_.del_vars(null) ;
        beta_.shed_rows(null) ;
        if (verbose) null.t().print("\tremoving variables") ;
      }
      
      // OPTIMALITY TESTING
      grad_ = - data_.XTy_ + set_.XTXA_ * beta_ ;
      // dual norm of gradient for inactive variables
      grd_norm = penalty_.elt_norm(grad_) - lambda_ ;
      // gradient for active variables
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
    if (current_it >= maxiter)  { status.back() = 1 ; }
    if (set_.size()  > maxfeat) { status.back() = 2 ; }
    if (!success)               { status.back() = 3 ; }
    
    // Preparing next value of the penalty
    if (status.back() >= 2) {
      break;
    } else {
      nzeros_ = join_cols(nzeros_, beta_/(data_.norm_X_(set_.A_) % lambda_factor_(set_.A_)));
      const_.push_back(data_.y_bar_ - as_scalar(dot(beta_, data_.X_bar_(set_.A_))));
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
  double loss = get_loss(), old_J = J_, old_D = D_ ;
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

