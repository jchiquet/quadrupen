/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#ifndef _GroupLasso_H
#define _GroupLasso_H

#include "GenericRegularizer.h"
#include "OptimizerL1.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

template <typename matrix>
class GroupLassoL1L2 : 
  public GenericRegularizer<matrix,Norm::L1L2>{
  public:
    
    using GenericRegularizer<matrix,Norm::L1L2>::intercept_ ;
    using GenericRegularizer<matrix,Norm::L1L2>::lambdas_   ;
    using GenericRegularizer<matrix,Norm::L1L2>::penalty_   ;
    using GenericRegularizer<matrix,Norm::L1L2>::data_ ;
    using GenericRegularizer<matrix,Norm::L1L2>::df_   ;
    using GenericRegularizer<matrix,Norm::L1L2>::lambda_factor_ ;
    using GenericRegularizer<matrix,Norm::L1L2>::get_lambda_seq ;
    
    GroupLassoL1L2(const RegressionData<matrix>&, const uvec&, const List&, const List&);
  
    double get_lambda_max() {return(penalty_.dual_norm(data_.XTy_, grp_sizes_));}
  
    List solution_path(const List&);
    
    uword nb_group() const { return grp_sizes_.n_elem ; }
    
    const sp_mat coefficients() const {
      return sp_mat(join_cols(iA_, jA_), nzeros_, data_.p_, lambdas_.size()) ; 
    }
    
    const sp_mat debiased_coefficients() const { 
      return sp_mat(join_cols(iA_, jA_), debiased_, data_.p_, lambdas_.size()) ; 
    }
    
    const vector<double>& intercept_debiased() const { return intercept_debiased_ ; }
    
    const sp_mat active_var() const { 
      return sp_mat(join_cols(iA_, jA_), vec(iA_.n_elem, fill::ones), data_.p_, lambdas_.size()) ; 
    }
    
    // Specific to Elastic-Net regularization
    GenericOptimizer<matrix,Norm::L1L2> solver_ ; // Solvers for L1 penalty
    ActiveSetGroup<matrix> set_ ; // Active set of variable and data
    double gamma_   ; // overall amount of l2 penalty
    vector<uvec> group_     ; // vector of current parameters
    uvec grp_sizes_ ; // vector of current parameters
    vec beta_       ; // vector of current parameters
    vec grad_       ; // vector of current gradient (smooth part)
    vec   nzeros_   ; // contains non-zero value of beta
    vec   debiased_ ; // contains the debiased non-zero value of beta
    vector<double >intercept_debiased_ ; // contains the debiased vector of intercept
    urowvec iA_     ; // contains row indices of the non-zero values
    urowvec jA_     ; // contains column indices of the non-zero values
    
    // Compute degrees of freedom for the current estimate
    double get_df() ;
    
};

template <typename matrix>
GroupLassoL1L2<matrix>::GroupLassoL1L2(
  const RegressionData<matrix>& data, const uvec& group_ind, const List& regParam, const List& control) :
  GenericRegularizer<matrix,Norm::L1L2>::GenericRegularizer(data, regParam) {

    // Vector of group and group sizes
    uvec grp = unique(group_ind) ;
    uword nb_grp =  grp.n_elem ;
    grp_sizes_.zeros(nb_grp-1) ;
    for (auto it = group_ind.begin(); it != group_ind.end(); it++) {
      grp_sizes_[*it - 1]++;
    }
    uword current_elt = 0 ;
    for (auto it = grp_sizes_.begin(); it != grp_sizes_.end(); it++) {
      group_.push_back(regspace<uvec>(current_elt, current_elt + *it - 1)) ;
      current_elt = current_elt + *it ;
    }
    
    // set the penalty to l1/l2
    penalty_ = Penalty<Norm::L1L2>() ;
    lambda_factor_ = as<vec>(regParam["lambda_factor"]) ;
    get_lambda_seq(get_lambda_max(), regParam) ;

    // Set up the optimizer
    solver_ = GenericOptimizer<matrix,Norm::L1L2>(penalty_);
    
    // Scale the structuring matrix according to main penalty factor and the amount of l2 penalty 
    gamma_   = as<double>(regParam["gamma"]) ;
    data_.scale_struct(sqrt(gamma_)*pow(lambda_factor_,-1)) ;
    data_.scale_regressors(lambda_factor_) ;
    
    // Initialize the active set, beta_ and gradient with starting coefficient
    vec beta0 = control["beta0"] ;
    // uvec A0 = find(beta0) ;
    // if (A0.is_empty()) {
    set_ = ActiveSetGroup(data_, grp_sizes_, as<bool>(control["usechol"])) ;
    grad_ = - data_.XTy_ ;
    // } else {
    //   set_  = ActiveSet(data_, A0, as<bool>(control["usechol"])) ;
    //   beta_ = beta0(A0) ;
    //   grad_ = - data_.XTy_ + set_.XTXA_ * beta_  ;
    // }
  
  }

template <typename matrix>
double GroupLassoL1L2<matrix>::get_df() {

  // TODO not correct at the moment

  double df = set_.size_grp() + data_.centered_ ;
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
List GroupLassoL1L2<matrix>::solution_path(const List& control) {
  
  // Parameters controlling the optimization
  const bool verbose(control["verbose"])      ; // verbosity level
  const double accuracy(control["threshold"]) ; // precision required
  const uword maxiter(control["maxiter"])     ; // max # of passes in the active set
  const uword maxfeat(control["maxfeat"])     ; // max # of variables activated
  const uword monitoring(control["monitor"])  ; // optimality monitor (0=none; 1=Grandvalet; 2=Fenchel)
  
  SolverType algorithm = FISTA; // Optimizer (default to FISTA)
  // if (as<std::string>(control["method"]) == "FISTA") algorithm = FISTA;

  // Variables monitoring the algorithm
  vector<double> gap, timing ; // timings and optimality measures
  vector<uword> status, iactive, ioptim    ; // convergence and # of inner/outer iterates

  // LAMBDA LOOP
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda_ : lambdas_) {
    if (verbose) {
      Rprintf("\n lambda_l1l2 = %f",lambda_) ;
      Rprintf("\n nb active groups = %i\n", set_.size_grp()) ;
    }
    
    // OPTIMIZER LOOP (FIX-LAMBDA VALUE): IDENTIFY THE ACTIVE SET AND SOLVE

    // variable associated with the highest violation of KKT conditions 
    vec optimality = penalty_.elt_norm(grad_, grp_sizes_) - lambda_;
    uword grp_in = optimality.index_max() ;
    double current_gap = std::max(0.0, optimality(grp_in)) ;
    uvec zeroed ;
    
    uword current_it = 0 ; bool success = true ; 
    while ((current_gap > accuracy) && (current_it <= maxiter)) {
      R_CheckUserInterrupt();
      current_it++;

      // ________________________________________________________________________
      // VARIABLE ACTIVATION IF APPLICABLE
      //
      if (set_.is_grp_in_[grp_in] == 0) { // Is var_in already in the active set?
        set_.add_group(grp_in, group_, data_) ;
        beta_.resize(beta_.size()+grp_sizes_(grp_in)) ; // update the vector of active parameters
        beta_.tail(grp_sizes_(grp_in)) = - 1e-3 * sign(grad_(group_[grp_in])) ;
        if (verbose) {Rprintf("\tnewly added group %i\n",grp_in);}
      } else if (verbose) {Rprintf("\talready in %i\n",grp_in);}

      // ________________________________________________________________________
      // OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      //
      // if (algorithm == FISTA) {
        ioptim.push_back(
          solver_.fista_groupwise(beta_, lambda_, data_, set_, 1e-10, 10000)
        );
      // } else { // QUADRA solver
        // try {
        //   ioptim.push_back(
        //     solver_.quadratic_enet(beta_, lambda_, data_, set_, 1e-5, 10000)
        //   );
        // } catch (std::runtime_error& error) {
        //   if (verbose > 0) {
        //     Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
        //   }
        //   success = false ;
        // }
      // }
      
      // OPTIMALITY TESTING
      optimality(set_.G_) = penalty_.elt_norm(grad_(set_.A_), grp_sizes_(set_.G_)) - lambda_ ;
      grp_in = optimality.index_max() ;
      current_gap = std::max(0.0, optimality(grp_in)) ;
      
    }
    if (verbose) Rprintf("\tcurrent gap = %f\n",current_gap) ;

    // Checking convergence status
    gap.push_back(current_gap) ;
    iactive.push_back(current_it) ;
    status.push_back(0) ;
    if (current_it >= maxiter) { status.back() = 1 ; }
    if (set_.size() > maxfeat) { status.back() = 2 ; }
    if (!success)              { status.back() = 3 ; }
    
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
      iA_ = join_rows(iA_, set_.A_.t()) ;
      jA_ = join_rows(jA_, df_.size()*ones<urowvec>(set_.size()) );
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

