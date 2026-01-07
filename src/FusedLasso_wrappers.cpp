#include <Rcpp.h>
using namespace Rcpp;

#include "FusedLasso_class.h"
#include "FusedLasso_data_struct.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List FusedLasso_cpp(
    const Environment &dataModel, // data structure
    const List  &tuningParam    , // List of tuning parameters
    const List  &control          // config of the optimisation 
  ) {

  const SEXP &R_XMat           = dataModel["X"]  ; // design matrix
  std::vector<double> y        = dataModel["y"]  ; // response vector
  List R_graph                 = dataModel["S"]  ; // Structuring matrix
  std::vector<double> penscale = dataModel["wx"] ;  // penalty weights
  std::vector<double> wObs     = dataModel["wy"] ; // responses to predictors vector
  // const arma::vec &normx    = dataModel["norm_X"] ; // norm of the predictors

  std::vector<double> lambda1Vec = tuningParam["l1"]   ; // vector of L1 penalties
  std::vector<double> lambda2Vec = tuningParam["l2"]   ; // scalar for the amount L2 penalty
  
  std::vector<double>  beta           = control["beta0"]         ;
  double mu0                          = control["mu0"]           ;
  const bool intercept                = control["intercept"]     ; 
  std::string penalty                 = control["pen_fused"]     ; 
  const unsigned int maxIterInner     = control["maxiterin"]     ;
  const unsigned int maxIterOuter     = control["maxiterout"]    ;
  const double accuracy               = control["accuracy"]      ;
  const unsigned int maxActivateVars  = control["maxactivation"] ;
  const unsigned int maxNonZero       = control["maxfeat"]       ;
  const bool adjustWithMax            = control["adjust"]        ;
  const int maxFusionLevel            = control["fusioncheck"]   ;
  const bool verbose                  = control["verbose"]       ;
  
  // Regression type (default to Gaussian)
  regEnum regType = GAUSSIAN;
  // if (family == "binomial") {regType = BINOMIAL;}
  
  // Penalty type (default to L1)
  penEnum penType = L1;
  if (penalty == "Huber") {penType = Huber;}
  if (penalty == "L2")    {penType = L2;}
  
  // Create the sparse matrix for X
  SparseMatrix X(R_XMat);

  // Import the connectivity graph
  Graph graph = Graph(R_graph["conn"], R_graph["weight"]);
  
  // Handle the intercept term
  if (intercept) {
    std::vector<double> ones(X.n, 1.0); 
    X.addColumn(ones) ;
    graph.addNode();
    penscale.push_back(1e-6) ;
    beta.push_back(mu0) ;
  }

  // Instantiate the main Fused-Lasso object
  FusedLasso fl(X, y, wObs, beta, penscale,  graph,
                maxIterInner, maxIterOuter, accuracy, 
                maxActivateVars, 0, 0, regType);

  // adjust the lambda1 for maximum value if necessary
  if (adjustWithMax) {
    list<int> exemptVars;
    for(unsigned int i = 0; i < penscale.size(); ++i) {
      if(penscale[i] < 1e-4) {
        exemptVars.push_back(i);   
      }
    }   
    double maxLambda1 = fl.findMaxLambda1(exemptVars);
    for(unsigned int i = 0; i < lambda1Vec.size(); ++i) {
      lambda1Vec[i] *= maxLambda1;
    }
  }

  // Vectors for monitoring optimization
  vector<bool> success;
  vector<int> outerIterNum;
  vector<int> innerIterNum;
  
  // Now, run
  SparseMatrix res = 
    fl.runAlgorithm(
      penType, maxFusionLevel, 
      lambda1Vec, lambda2Vec, maxNonZero, 
      success, outerIterNum, innerIterNum, verbose
  );

  return List::create(
    Named("tuning_param") = List::create(
      Named("l1") = lambda1Vec,
      Named("l2") = lambda2Vec 
    ),
    Named("beta")        = wrap(res.todgCMatrix()),
    Named("monitoring")  = 
      List::create(
        Named("it_active")      = outerIterNum,
        Named("it_optim")       = innerIterNum ,
        Named("convergence")    = success
      )
  );

}
