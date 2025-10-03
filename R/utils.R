## ======================================================
## GENERATE A GRID OF PENALTY IF NONE HAS BEEN PROVIDED
get.lambda1.l1 <- function(xty,nlambda1,min.ratio) {
  lmax <- max(abs(xty))
  return(10^seq(log10(lmax), log10(min.ratio*lmax), len=nlambda1))
}

get.lambda1.li <- function(xty,nlambda1,min.ratio) {
  lmax <- sum(abs(xty))
  return(10^seq(log10(lmax), log10(min.ratio*lmax), len=nlambda1))
}

## ======================================================
## RECOVER THE LIST OF DEFAULT OPTIONAL ARGUMENTS
default.args <- function(penalty,n,p,user) {
  switch(penalty,
         "elastic.net" = list(
           lambda1   = NULL,
           lambda2   = 0.01,
           penscale  = rep(1,p),
           struct    = NULL,
           intercept = TRUE,
           normalize = TRUE,
           naive     = FALSE,
           nlambda1  = ifelse(is.null(user$lambda1),100,length(user$lambda1)),
           min.ratio = ifelse(n<p,0.01,5e-3),
           max.feat  = min(4*n,p),
           beta0     = NULL,
           control   = list(),
           checkargs = TRUE),

         "lasso"     = list(
           lambda1   = NULL,
           penscale  = rep(1,p),
           intercept = TRUE,
           normalize = TRUE,
           nlambda1  = ifelse(is.null(user$lambda1),100,length(user$lambda1)),
           min.ratio = ifelse(n<p,0.01,5e-3),
           max.feat  = min(n,p),
           beta0     = NULL,
           control   = list(),
           checkargs = TRUE),
         
         "ridge" = list(
           lambda2    = NULL,
           struct     = NULL,
           intercept  = TRUE,
           normalize  = TRUE,
           nlambda2   = 100 ,
           lambda.min = ifelse(n<=p,0.01,1e-4),
           lambda.max = 100 ,
           control    = list(),
           checkargs  = TRUE),

         "bounded.reg" = list(
           lambda1   = NULL,
           lambda2   = 0.01,
           penscale  = rep(1,p),
           struct    = NULL,
           intercept = TRUE,
           normalize = TRUE,
           naive     = FALSE,
           nlambda1  = ifelse(is.null(user$lambda1),100,length(user$lambda1)),
           min.ratio = ifelse(n<p,0.01,5e-3),
           max.feat  = min(4*n,p),
           control   = list(),
           checkargs = TRUE)
         )
}

## ====================================================================
## STANDARDIZE THE PREDICTOR (NEEDED FOR CROSS-VALIDATION PURPOSES)
standardize <- function(x,y,intercept,normalize,penscale,zero=.Machine$double.eps) {

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
  ## Computing marginal correlation
  if (intercept) {
    xty   <- drop(crossprod(y-ybar,sweep(x,2L,xbar)))
  } else {
    xty   <- drop(crossprod(y,x))
  }

  return(list(xbar=xbar, ybar=ybar, normx=normx, normy=normy, xty=xty, x=x))
}

ctrl_default <- function(d)
  list(verbose    = 1, # default control options
      timer       =  FALSE,
      max.iter    = max(500, d),
      method      = "quadra",
      threshold   = 1e-7,
      monitor     = 0,
      bulletproof = TRUE,
      usechol     = TRUE
  )



