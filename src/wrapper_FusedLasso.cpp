// [[Rcpp::depends(RcppArmadillo)]]

// Include Armadillo / Rcpp / R to C/C++ basics
#include "RcppArmadillo.h"
#include "FusedLasso/FusedLasso_class.h"
#include "FusedLasso/FusedLasso_data_struct.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

vec get_lambda1(SEXP LAMBDA1, uword n_lambda, const double min_ratio, double lmax) {
  vec lambda1 ;
  if (LAMBDA1 != R_NilValue) {
    lambda1  = as<vec>(LAMBDA1)  ;
  } else {
    lambda1 = exp10(linspace(log10(lmax), log10(min_ratio*lmax), n_lambda)) ;
  }
  return(lambda1);
}

// [[Rcpp::export]]
Rcpp::List FusedLasso_cpp(
    const Environment &dataModel   , // data structure
    const bool        &intercept   , // boolean for intercept
    const List        &regParam    , // config of the optimisation 
    const List        &controlFit  // config of the optimisation   
) {
  const uword n                = dataModel["n"]  ; // sample size
  const uword p                = dataModel["d"]  ; // number of features
  const SEXP &R_XMat           = dataModel["X"]  ; // design matrix
  std::vector<double> y        = dataModel["y"]  ; // response vector
  List R_graph                 = dataModel["S"]  ; // Structuring matrix
  std::vector<double> wObs     = dataModel["wy"] ; // responses to predictors vector

  const SEXP R_LAMBDA1         = regParam["l1"]         ; // vector of L1 penalties ; if NULL, automatically set
  std::vector<double> wlambda1 = regParam["l1_weights"] ; // l1-penalty weights
  double lambda2               = regParam["l2"]   ; // scalar for the amount L2 penalty
  uword n_lambda               = regParam["n_lambda1"]  ; // # of l1-penalty levels
  const double min_ratio       = regParam["min_ratio"]  ; // minimum penlaty value as a ratio of lambda1 max

  std::vector<double>  beta0          = controlFit["beta0"]         ;
  double mu0                          = controlFit["mu0"]           ;
  std::string penalty                 = controlFit["pen_fused"]     ; 
  const unsigned int maxIterInner     = controlFit["maxiterin"]     ;
  const unsigned int maxIterOuter     = controlFit["maxiterout"]    ;
  const double accuracy               = controlFit["accuracy"]      ;
  const unsigned int maxActivateVars  = controlFit["maxactivation"] ;
  const unsigned int maxNonZero       = controlFit["maxfeat"]       ;
  const int maxFusionLevel            = controlFit["fusioncheck"]   ;
  const bool verbose                  = controlFit["verbose"]       ;
  const bool normalize                = controlFit["normalize"]       ;
  
  // Regression type (default to Gaussian)
  regEnum regType = GAUSSIAN;

  // Penalty type (default to L1)
  penEnum penType = L1;
  if (penalty == "Huber") {penType = Huber;}
  if (penalty == "L2")    {penType = L2;}

  // Create the sparse data matrix for X
  arma::sp_mat sp_x(as<sp_mat>(R_XMat)) ;
  vec normx ;
  if (normalize) {
    normx = sqrt(trans(sum(square(sp_x),0)));
    for (uword i=0; i<p; i++) {
      sp_x.col(i) /= normx(i);
    }
  } else {
    normx = ones(p);
  }
  SparseMatrix X(sp_x);

  // Import the connectivity graph
  Graph graph = Graph(R_graph["conn"], R_graph["weight"]);
  
  // Handle the intercept term
  if (intercept) {
    vector<double> ones(n, 1.0); 
    X.addColumn(ones) ;
    graph.addNode();
    wlambda1.push_back(1e-6) ;
    beta0.push_back(mu0) ;
  }
  
  // Instantiate the main Fused-Lasso object
  FusedLasso fl(X, y, wObs, beta0, wlambda1,  graph,
                maxIterInner, maxIterOuter, accuracy, 
                maxActivateVars, 0, 0, regType);

  // Vectors of penalties
  list<int> exemptVars;
  for(unsigned int i = 0; i < wlambda1.size(); ++i) {
    if(wlambda1[i] < 1e-4) {
      exemptVars.push_back(i);   
    }
  }
  double maxLambda1 = fl.findMaxLambda1(exemptVars);
  vec lambda_l1 = get_lambda1(R_LAMBDA1, n_lambda, min_ratio, maxLambda1) ;
  n_lambda = lambda_l1.n_elem ; // # of penalty levels
    vector<double> lambda1Vec = arma::conv_to < std::vector<double> >::from(lambda_l1) ;
  vector<double> lambda2Vec(lambda1Vec.size(), lambda2) ;

  // OPTIMIZATION BY COORDINATE DESCENT
  vector<bool> success;      // Vectors for monitoring optimization
  vector<int> outerIterNum;
  vector<int> innerIterNum;
  SparseMatrix res = 
    fl.runAlgorithm(
      penType, maxFusionLevel, 
      lambda1Vec, lambda2Vec, maxNonZero, 
      success, outerIterNum, innerIterNum, verbose
  );
  
  // Preparing R output
  arma::sp_mat beta = Rcpp::as<arma::sp_mat>(res.todgCMatrix()) ;
  beta = beta.t() ;
  vec mu = arma::zeros(lambda1Vec.size()) ;
  if (intercept) {
    mu = beta.col(p) ;
    beta = beta.cols(0, p-1) ;
  }
  
  // degrees of freedom : number of non-zero values
  arma::vec df(beta.n_rows) ;
  for (uword i=0; i<beta.n_rows; i++) {
    rowvec row = beta.row(i).as_dense() ;
    vec val = nonzeros(unique(row.t())) ;
    df(i) = val.n_elem ;
  }

  // Normalizing beta back to original scale
  beta = beta * arma::diagmat(1/normx) ;

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = NumericVector(lambda1Vec.begin(), lambda1Vec.begin()+beta.n_rows),
      Named("l2") = lambda2
    ),
    Named("beta")       = beta,
    Named("active")     = spones(beta),
    Named("mu")         = mu,
    Named("df")         = df,
    Named("monitoring") = 
      List::create(
        Named("it_active")      = outerIterNum,
        Named("it_optim")       = innerIterNum ,
        Named("convergence")    = success
      )
  );

}
