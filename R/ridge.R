#' Fit a linear model with a structured ridge regularization
#'
#' Adjust a linear model with ridge regularization (possibly
#' structured \eqn{\ell_2}{l2}-norm). The solution path is computed
#' at a grid of values for the \eqn{\ell_2}{l2}-penalty. See details
#' for the criterion optimized.
#'
#' @param x matrix of features. Do NOT include intercept. When
#' normalized os \code{TRUE}, coefficients will then be rescaled to
#' the original scale.
#'
#' @param y response vector.
#'
#' @param lambda2 sequence of decreasing \eqn{\ell_2}{l2}-penalty
#' levels. If \code{NULL} (the default), a vector is generated with
#' \code{nlambda2} entries.
#'
#' @param struct matrix structuring the coefficients, possibly
#' sparsely encoded. Must be at least positive semi-definite (this is
#' checked internally if the \code{checkarg} argument is
#' \code{TRUE}). If \code{NULL} (the default), the identity matrix is
#' used. See details below.
#'
#' @param intercept logical; indicates if an intercept should be
#' included in the model. Default is \code{TRUE}.
#'
#' @param normalize logical; indicates if variables should be
#' normalized to have unit L2 norm before fitting.  Default is
#' \code{TRUE}.
#'
#' @param nlambda2 integer that indicates the number of values to put
#' in the \code{lambda2} vector.  Ignored if \code{lambda2} is
#' provided.
#'
#' @param lambda.min the minimal amount of penalty used to generated
#' the vector \code{lambda2}. Ignored if \code{lambda2} is provided.
#'
#' @param lambda.max the maximal amount of penalty used to generated
#' the vector \code{lambda2}. Ignored if \code{lambda2} is provided.
#'
#' @param control list of argument controlling low level options of
#' the algorithm --use with care and at your own risk-- :
#' \itemize{%
#'
#' \item{\code{verbose}: }{integer; activate verbose mode --this one
#' is not too much risky!-- set to \code{0} for no output; \code{1}
#' for warnings only, and \code{2} for tracing the whole
#' progression. Default is \code{1}. Automatically set to \code{0}
#' when the method is embedded within cross-validation or stability
#' selection.}
#'
#' \item{\code{timer}: }{logical; use to record the timing of the
#' algorithm. Default is \code{FALSE}.}
#'
#' }
#'
#' @param checkargs logical; should arguments be checked to
#' (hopefully) avoid internal crashes? Default is
#' \code{TRUE}. Automatically set to \code{FALSE} when calls are made
#' from cross-validation or stability selection procedures.
#'
#' @return an object with class \code{quadrupen}, see the
#' documentation page \code{\linkS4class{quadrupen}} for details.
#'
#' @note The optimized criterion is the following: \if{latex}{\deqn{%
#' \hat{\beta}_{\lambda_2} = \arg \min_{\beta} \frac{1}{2} (y - X
#' \beta)^T (y - X \beta) + \frac{\lambda_2}{2} \beta^T S \beta, }}
#' \if{html}{\out{ <center> &beta;<sup>hat</sup>
#' <sub>&lambda;<sub>2</sub></sub> = argmin<sub>&beta;</sub> 1/2
#' RSS(&beta) + + &lambda;/2 <sub>2</sub> &beta;<sup>T</sup> S
#' &beta;, </center> }} \if{text}{\deqn{beta.hat(lambda2) =
#' argmin_beta 1/2 RSS(beta) + lambda2 beta' S beta,}} where the
#' \eqn{\ell_2}{l2} structuring positive semidefinite matrix
#' \eqn{S}{S} is provided via the \code{struct} argument (possibly of
#' class \code{Matrix}).
#'
#' @seealso See also \code{\linkS4class{quadrupen}},
#' \code{\link{plot.quadrupen}} and \code{\link{crossval}}.
#' @name ridge
#' @rdname ridge
#' @keywords models, regression
#' 
#' @examples
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
#' cor <- 0.75
#' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
#' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
#' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
#' diag(Sigma) <- 1
#' n <- 50
#' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
#' y <- 10 + x %*% beta + rnorm(n,0,10)
#'
#'
#' beta.lasso <- slot(crossval(x,y, penalty="lasso", mc.cores=2) , "beta.min")
#'
#' cat("\nFalse positives for the Lasso:", sum(sign(beta) != sign(beta.lasso)))
#' cat("\nDONE.\n")
#'
#'
#' @export
ridge <- function(x,
                  y,
                  lambda2    = NULL,
                  struct     = NULL,
                  intercept  = TRUE,
                  normalize  = TRUE,
                  nlambda2   = 100 ,
                  lambda.min = ifelse(n<=p,1e-2,1e-4),
                  lambda.max = 100,
                  control    = list(),
                  checkargs  = TRUE) {

  p <- ncol(x) # problem size
  n <- nrow(x) # sample size

  ## ===================================================
  ## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
  if(!inherits(x, c("matrix")))
      x <- as.matrix(x)
  if (checkargs) {
    if(any(is.na(x)))
      stop("NA value in x not allowed.")
    if(!is.numeric(y))
      stop("y has to be of type 'numeric'")
    if(n != length(y))
      stop("x and y have not correct dimensions")
    if (!is.null(lambda2)) {
      if(any(lambda2 <= 0) | lambda.min <=0)
        stop("entries in lambda2 must all be postive.")
    }
    if (!is.null(struct)) {
      if (ncol(struct) != p | ncol(struct) != p)
          stop("struct must be a (square) positive semidefinite matrix.")
      if (any(eigen(struct,only.values=TRUE)$values<0))
          stop("struct must be a (square) positive semidefinite matrix.")
      if(!inherits(struct, c("dgCMatrix", "matrix")))
          struct <- as(struct, "dgCMatrix")
    }
  }

  ## ============================================
  ## RECOVERING LOW LEVEL OPTIONS
  quadra <- TRUE
  if (!is.null(control$method)) {
    if (control$method != "quadra") {
      quadra <- FALSE
    }
  }
  ctrl <- list(verbose      = 1, # default control options
               timer        =  FALSE)
  ctrl[names(control)] <- control # overwritten by user specifications
  if (ctrl$timer) {r.start <- proc.time()}

  ## Cholesky decomposition of the (possibly sparsely encoded) structuring matrix
  if (is.null(struct)) {
    C <- diag(rep(1,p))
  } else {
    C <- as.matrix(chol(struct))
  }

  ## ======================================================
  ## STARTING C++ CALL TO ENET_LS
  if (ctrl$timer) {cpp.start <- proc.time()}
  out <- ridge_cpp(
    x            ,
    y            ,
    C            ,
    lambda2      ,
    nlambda2     ,
    lambda.min   ,
    lambda.max   ,
    intercept    ,
    normalize    ,
    rep(1,n)     ,
    ctrl$verbose)
  coefficients <- Matrix(out$coefficients)
  ## END OF CALL
  if (ctrl$timer) {
    internal.timer <- (proc.time() - cpp.start)[3]
    external.timer <- (proc.time() - r.start)[3]
  } else {
    internal.timer <- NULL
    external.timer <- NULL
  }

  ## ======================================================
  ## BUILDING THE QUADRUPEN OBJECT
  monitoring  <- list(external.timer = external.timer   ,
                      internal.timer = internal.timer   )
  dimnames(coefficients)[[1]] <- round(c(out$lambda2),3)
  if (is.null(colnames(x))) {
    dimnames(coefficients)[[2]] <- 1:p
  } else {
    dimnames(coefficients)[[2]] <- colnames(x)
  }
  mu <- drop(out$mu)
  df <- drop(out$df)

  ## FITTED VALUES AND RESIDUALS...
  if (intercept) {
    fitted <- sweep(tcrossprod(x,coefficients),2L,-mu,check.margin=FALSE)
    df <- df + 1
  } else {
    mu <- 0
    fitted <- tcrossprod(x,coefficients)
  }
  residuals <- apply(fitted,2,function(y.hat) y - y.hat)
  r.squared <- 1 - colSums(residuals^2)/ifelse(intercept,sum((y - mean(y))^2),sum(y^2))

  return(new("quadrupen",
             coefficients = coefficients   ,
             active.set   = Matrix(1,nrow=p,ncol=p),
             intercept    = intercept      ,
             mu           = mu             ,
             normx        = drop(out$normx),
             fitted       = fitted         ,
             residuals    = residuals      ,
             df           = df             ,
             r.squared    = r.squared      ,
             penscale     = rep(1,p)       ,
             penalty      = "ridge"        ,
             naive        = NULL           ,
             lambda1      = 0              ,
             lambda2      = c(out$lambda2) ,
             monitoring   = monitoring     ,
             control      = ctrl))

}
