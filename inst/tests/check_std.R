library(inline)
library(RcppArmadillo)
library(testthat)
library(Matrix)

standardize_R <- function(x,y,intercept,normalize,penscale,zero=.Machine$double.eps) {

  n <- length(y)
  p <- ncol(x)
  ## ============================================
  ## INTERCEPT AND NORMALIZATION TREATMENT
  if (intercept) {
    xbar <- colMeans(x)
    ybar <- mean(y)
  } else {
    xbar <- rep(0,p)
    ybar <- 0
  }

  ## ============================================
  ## NORMALIZATION
  if (normalize) {
    normx <- sqrt(drop(colSums(x^2)- n*xbar^2))
    if (any(normx < zero)) {
      warning("A predictor has no signal: you should remove it.")
      normx[abs(normx) < zero] <- 1 ## dirty way to handle 0/0
    }
    ## normalizing the predictors...
    x <- sweep(x, 2L, normx, "/", check.margin = FALSE)
    ## xbar is scaled to handle internaly the centering of X for
    ## sparsity purpose
    xbar <- xbar/normx
  } else {
    normx <- rep(1,p)
  }
  normy <- sqrt(sum(y^2))

  ## and now normalize predictors according to penscale value
  if (any(penscale != 1)) {
    x <- sweep(x, 2L, penscale, "/", check.margin=FALSE)
    xbar <- xbar/penscale
  }
  ## Computuping marginal correlation
  if (intercept) {
    xty   <- drop(crossprod(y-ybar,sweep(x,2L,xbar)))
  } else {
    xty   <- drop(crossprod(y,x))
  }

  return(list(xbar=xbar, ybar=ybar, normx=normx, normy=normy, xty=xty, x=x))
}

code_dense <- '
  using namespace Rcpp;
  using namespace arma;

  mat  x         = as<mat> (X)         ;
  vec  y         = as<vec> (Y)         ;
  bool intercept = as<bool>(INTERCEPT) ;
  bool normalize = as<bool>(NORMALIZE) ;
  vec  penscale  = as<vec> (PENSCALE)  ;

  rowvec xbar     ;
  double ybar  ;
  rowvec normx    ;
  double normy ;
  vec xty      ;

  uword n = x.n_rows;
  uword p = x.n_cols;

  if (intercept == 1) {
    xbar = mean(x, 0);
    ybar = mean(y) ;
  } else {
    xbar = zeros(1,p) ;
    ybar = 0;
  }

  if (normalize == 1) {
    normx = sqrt(sum(square(x),0) - n * square(xbar));
    for (int i=0; i<p; i++) {
      x.col(i) = x.col(i) / normx(i);
    }
    xbar /= normx ;
  } else {
    normx = ones(1, p);
  }
  normy = sqrt(sum(square(y))) ;

  if (any(penscale != 1)) {
    for (int i=0; i<p; i++) {
       x.col(i) = x.col(i)/ penscale(i) ;
    }
    xbar /= penscale;
  }

  if (intercept == 1) {
    xty = trans(trans(y-ybar) * x) ;
    for (int i=0;i<p;i++) {
       xty(i) -=  sum(y-ybar) * xbar(i);
    }
  } else {
    xty = trans(y.t()*x) ;
  }

  return List::create(Named("xbar")  = xbar ,
		      Named("ybar")  = ybar ,
		      Named("normx") = normx,
		      Named("normy") = normy,
		      Named("xty")   = xty  ,
		      Named("x")     = wrap(x)   ) ;
'

standardize_cpp <- cxxfunction(signature(X="matrix",Y="numeric",INTERCEPT="boolean",
                                         NORMALIZE="boolean",PENSCALE="numeric"),
                               plugin="RcppArmadillo", body=code_dense)


x <- matrix(rnorm(2000*1000),2000,1000)
y <- rnorm(2000)
weights <- rep(1,ncol(x))
print(system.time(dense.R   <- standardize_R(x,y,TRUE,TRUE,weights)))
print(system.time(dense.cpp <- standardize_cpp(x,y,TRUE,TRUE,weights)))

dense.cpp$xbar <- as.numeric(dense.cpp$xbar)
dense.cpp$normx <- as.numeric(dense.cpp$normx)
dense.cpp$xty <- as.numeric(dense.cpp$xty)

expect_that(dense.cpp, is_equivalent_to(dense.R))

x[sample(1:length(x),2000*800)] <- 0
x <- Matrix(x, sparse=TRUE)

code_sparse <- '
  using namespace Rcpp;
  using namespace arma;

  sp_mat  x      = as<sp_mat> (X)         ;
  vec  y         = as<vec> (Y)         ;
  bool intercept = as<bool>(INTERCEPT) ;
  bool normalize = as<bool>(NORMALIZE) ;
  vec  penscale  = as<vec> (PENSCALE)  ;

  rowvec xbar  ;
  double ybar  ;
  rowvec normx ;
  double normy ;
  vec xty      ;

  uword n = x.n_rows;
  uword p = x.n_cols;

  if (intercept == 1) {
    xbar = rowvec(mean(x, 0));
    ybar = mean(y) ;
  } else {
    xbar = zeros(1,p) ;
    ybar = 0;
  }

  if (normalize == 1) {
    normx = sqrt(sum(square(x),0) - n * square(xbar));
    for (int i=0; i<p; i++) {
      x.col(i) /= normx(i);
    }
    xbar /= normx ;
  } else {
    normx = ones(1, p);
  }
  normy = sqrt(sum(square(y))) ;

  if (any(penscale != 1)) {
    for (int i=0; i<n; i++) {
       x.row(i) /= penscale ;
    }
    xbar /= penscale;
  }

  if (intercept == 1) {
    xty = trans(trans(y-ybar) * x) ;
    for (int i=0;i<p;i++) {
       xty(i) -=  sum(y-ybar) * xbar(i);
    }
  } else {
    xty = trans(y.t()*x) ;
  }

  return List::create(Named("xbar")  = xbar ,
		      Named("ybar")  = ybar ,
		      Named("normx") = normx,
		      Named("normy") = normy,
		      Named("xty")   = xty  ,
		      Named("x")     = wrap(x)) ;
'

standardize_cpp_sparse <- cxxfunction(signature(X="dgCMatrix",Y="numeric",INTERCEPT="boolean",
                                                NORMALIZE="boolean",PENSCALE="numeric"),
                                      plugin="RcppArmadillo",
                                      include= "#include <RcppArmadilloExtensions/spmat.h>",
                                      body=code_sparse)

print(system.time(sparse.R   <- standardize_R(x,y,TRUE,TRUE,rep(1,1000))))
print(system.time(sparse.cpp <- standardize_cpp_sparse(x,y,TRUE,TRUE,rep(1,1000))))

sparse.cpp$xbar <- as.numeric(sparse.cpp$xbar)
sparse.cpp$normx <- as.numeric(sparse.cpp$normx)
sparse.cpp$xty <- as.numeric(sparse.cpp$xty)
sparse.cpp$x <- as.matrix(sparse.cpp$x)
sparse.R$x <- as(sparse.R$x, "matrix")

expect_that(sparse.cpp, is_equivalent_to(sparse.R))
