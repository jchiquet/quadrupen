##' Fit a linear model with group-lasso (either l1/l2 or l1/l-inf) regularization
##'
##' Adjust a linear model with group-lasso regularization, that is a
##' mixture of either a (possibly weighted)
##' \eqn{\ell_1/\ell_2}{l1/l2}- or
##' \eqn{\ell_1/\ell_\infty}{l1/linf}-norm, and a (possibly
##' structured) \eqn{\ell_2}{l2}-norm (ridge-like). The solution path
##' is computed at a grid of values for the
##' \eqn{\ell_1/\ell_q}{l1/lq}-penalty. See details for the criterion
##' optimized.
##'
##' @param x matrix of features, possibly sparsely encoded
##' (experimental). Do NOT include intercept. When normalized os
##' \code{TRUE}, coefficients will then be rescaled to the original
##' scale.
##'
##' @param y response vector.
##'
##' @param group vector of interers indicating group belonging. Must
##' match the number fo column in \code{x}. Must be SORTED integers
##' starting from 1.
##'
##' @param type string indicating wether the
##' \eqn{\ell_1/\ell_2}{l1/l2} or the
##' \eqn{\ell_1/\ell_\infty}{l1/linf} group-Lasso must be fitted.
##'
##' @param lambda1 sequence of decreasing \eqn{\ell_1}{l1}-penalty
##' levels. If \code{NULL} (the default), a vector is generated with
##' \code{nlambda1} entries, starting from a guessed level
##' \code{lambda1.max} where only the intercept is included, then
##' shrunken to \code{min.ratio*lambda1.max}.
##'
##' @param lambda2 real scalar; tunes the \eqn{\ell_2}{l2}
##' penalty. Default is 0.01. Set to 0 to recover the usual
##' group-Lasso.
##'
##' @param penscale vector with real positive values that weight the
##' \eqn{\ell_1}{l1}-penalty of each feature.
##'
##' @param struct matrix structuring the coefficients, possibly
##' sparsely encoded. Must be at least positive semidefinite (this is
##' checked internally if the \code{checkarg} argument is
##' \code{TRUE}). If \code{NULL} (the default), the identity matrix is
##' used. See details below.
##'
##' @param intercept logical; indicates if an intercept should be
##' included in the model. Default is \code{TRUE}.
##'
##' @param normalize logical; indicates if variables should be
##' normalized to have unit L2 norm before fitting.  Default is
##' \code{TRUE}.
##'
##' @param naive logical; Compute either 'naive' of classic version:
##' as defined for the Elastic-net in Zou and Hastie (2006), the
##' vector of parameters is rescaled by a coefficient
##' \code{(1+lambda2)} when \code{naive} equals \code{FALSE}.  No
##' rescaling otherwise.  Default is \code{FALSE}.
##'
##' @param nlambda1 integer that indicates the number of values to put
##' in the \code{lambda1} vector.  Ignored if \code{lambda1} is
##' provided.
##'
##' @param min.ratio minimal value of \eqn{\ell_1}{l1}-part of the
##' penalty that will be tried, as a fraction of the maximal
##' \code{lambda1} value. A too small value might lead to unstability
##' at the end of the solution path corresponding to small
##' \code{lambda1}.  The default value tries to avoid this, adapting
##' to the '\eqn{n<p}{n<p}' context. Ignored if \code{lambda1} is
##' provided.
##'
##' @param max.feat integer; limits the number of features ever to
##' enter the model; i.e., non-zero coefficients for the Lasso: the
##' algorithm stops if this number is exceeded and \code{lambda1} is
##' cutted at the corresponding level. Default is
##' \code{min(nrow(x),ncol(x))}. Use with care, as it considerably
##' changes the computation time.
##'
##' @param beta0 a starting point for the vector of parameter. The
##' default is a vector of zeros. May save time in some situation.
##'
##' @param control list of argument controlling low level options of
##' the algorithm --use with care and at your own risk-- :
##' \itemize{%
##'
##' \item{\code{verbose}: }{integer; activate verbose mode --this one
##' is not too much risky!-- set to \code{0} for no output; \code{1}
##' for warnings only, and \code{2} for tracing the whole
##' progression. Default is \code{1}. Automatically set to \code{0}
##' when the method is embedded within cross-validation or stability
##' selection.}
##'
##' \item{\code{timer}: }{logical; use to record the timing of the
##' algorithm. Default is \code{FALSE}.}
##'
##' \item{\code{max.iter}: }{the maximal number of iteration used to
##' solve the problem for a given value of lambda1. Default is 500.}
##'
##' \item{\code{method}: }{a string for the underlying solver
##' used. Either \code{"quadra"}, \code{"pathwise"} or
##' \code{"fista"}. Default is \code{"quadra"}.}
##'
##' \item{\code{threshold}: }{a threshold for convergence. The
##' algorithm stops when the optimality conditions are fulfill up to
##' this threshold. Default is \code{1e-7} for \code{"quadra"} and
##' \code{1e-2} for the first order methods.}
##'
##' \item{\code{monitor}: }{indicates if a monitoring of the
##' convergence should be recorded, by computing a lower bound between
##' the current solution and the optimum: when \code{'0'} (the
##' default), no monitoring is provided; when \code{'1'}, the bound
##' derived in Grandvalet et al. is computed; when \code{'>1'}, the
##' Fenchel duality gap is computed along the algorithm.}
##' }
##'
##' @param checkargs logical; should arguments be checked to
##' (hopefully) avoid internal crashes? Default is
##' \code{TRUE}. Automatically set to \code{FALSE} when calls are made
##' from cross-validation or stability selection procedures.
##'
##' @return an object with class \code{quadrupen}, see the
##' documentation page \code{\linkS4class{quadrupen}} for details.
##'
##' @note The optimized criterion is the following: \if{latex}{\deqn{%
##' \hat{\beta}_{\lambda_1} = \arg \min_{\beta} \frac{1}{2} (y - X
##' \beta)^T (y - X \beta) + \lambda_1 \|D \beta \|_{1}, }}
##' \if{html}{\out{ <center> &beta;<sup>hat</sup>
##' <sub>&lambda;<sub>1</sub></sub> =
##' argmin<sub>&beta;</sub> 1/2 RSS(&beta) + &lambda;<sub>1</sub>
##' &#124; D &beta; &#124;<sub>1</sub>, </center> }}
##' \if{text}{\deqn{beta.hat(lambda1) = argmin_beta 1/2
##' RSS(beta) + lambda1 |D beta|1,}} where
##' \eqn{D}{D} is a diagonal matrix, whose diagonal terms are provided
##' as a vector by the \code{penscale} argument.
##'
##' @seealso See also \code{\linkS4class{quadrupen}},
##' \code{\link{plot.quadrupen}} and \code{\link{crossval}}.
##' @name group.lasso
##' @rdname group.lasso
##' @keywords models, regression
##'
##'
##' @export
group.lasso <- function(x,
                        y,
                        group,
                        type        = 2, ## could be "inf" or 2, default is 2
                        lambda1   = NULL,
                        lambda2   = 0.01,
                        penscale  = switch(as.character(type), "2" = rep(sqrt(pk), pk), "inf" = rep(1,p)),
                        struct    = NULL,
                        intercept = TRUE,
                        normalize = TRUE,
                        naive     = FALSE,
                        nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                        min.ratio = ifelse(n<=p,1e-2,1e-4),
                        max.feat  = ifelse(lambda2<1e-2,min(n,p),min(4*n,p)),
                        beta0     = rep(0,p),
                        control   = list(),
                        checkargs = TRUE) {

  p <- ncol(x) # problem size
  n <- nrow(x) # sample size
  pk <- tabulate(group) # group size

  ## ===================================================
  ## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
  if (checkargs) {
    if (is.data.frame(x))
      x <- as.matrix(x)
    if(!inherits(x, c("matrix", "dgCMatrix")))
      stop("x has to be of class 'matrix' or 'dgCMatrix'.")
    if(any(is.na(x)))
      stop("NA value in x not allowed.")
    if(!is.numeric(y))
      stop("y has to be of type 'numeric'")
    if(n != length(y))
      stop("x and y have not correct dimensions")
    if(length(penscale) != p)
      stop("penscale must have ncol(x) entries")
    if (any(penscale <= 0))
      stop("weights in penscale must be positive")
    if(!inherits(lambda2, "numeric") | length(lambda2) > 1)
      stop("lambda2 must be a scalar.")
    if(lambda2 < 0)
      stop("lambda2 must be a non negative scalar.")
    if (!is.null(lambda1)) {
      if(any(lambda1 <= 0))
        stop("entries inlambda1 must all be postive.")
      if(is.unsorted(rev(lambda1)))
        stop("lambda1 values must be sorted in decreasing order.")
      if (length(lambda1)>1 & !is.null(beta0))
        warning("providing beta0 for a serie of l1 penalties mught be inefficient.")
    }
    if(min.ratio < 0)
        stop("min.ratio must be non negative.")
    if (!is.null(struct)) {
      if (ncol(struct) != p | ncol(struct) != p)
          stop("struct must be a (square) positive semidefinite matrix.")
      if (any(eigen(struct,only.values=TRUE)$values<0))
          stop("struct must be a (square) positive semidefinite matrix.")
      if(!inherits(struct, "dgCMatrix"))
          struct <- as(struct, "dgCMatrix")
    }
    ## Check the group vector
    if(length(group) != p)
        stop("group length does not match ncol(x)")
    if(is.unsorted(group))
        stop("the group vector must be sorted (do the manual permutation of the columns of x please)")
    if(any(diff(unique(group)) != 1))
        stop("the groups must be labeled from 1 to nb groups with increment of 1.")

    if (!is.null(beta0)) {
      beta0 <- as.numeric(beta0)
      if (length(beta0) != p)
        stop("beta0 must be a vector with p entries.")
    }
    if (length(max.feat)>1)
        stop("max.feat must be an integer.")
    if(is.numeric(max.feat) & !is.integer(max.feat))
        max.feat <- as.integer(max.feat)
  }

  ## ============================================
  ## RECOVERING LOW LEVEL OPTIONS
  quadra <- TRUE
  if (!is.null(control$method)) {
    if (control$method != "quadra") {
      quadra <- FALSE
    }
  }
  ctrl <- list(verbose   = FALSE, # default control options
               timer     = FALSE,
               zero      = .Machine$double.eps,
               max.iter  = 500,
               method    = ifelse(type=='inf', "quadra","fista"),
               threshold = ifelse(quadra, 1e-7, 1e-2))

  ctrl[names(control)] <- control # overwritten by user specifications
  if (ctrl$timer) {r.start <- proc.time()}

  ## ======================================================
  ## STARTING C++ CALL TO group_l1linf
  if (ctrl$timer) {cpp.start <- proc.time()}
  out <- .Call("group_lasso",
               beta0        ,
               x            ,
               y            ,
               switch(type, "inf" = 1, 2),
               pk           ,
               struct       ,
               lambda1      ,
               nlambda1     ,
               min.ratio    ,
               penscale     ,
               lambda2      ,
               intercept    ,
               normalize    ,
               rep(1,n)     ,
               naive        ,
               ctrl$thresh  ,
               ctrl$max.iter,
               max.feat     ,
               switch(ctrl$method,
                      quadra   = 0,
                      pathwise = 1,
                      fista    = 2, 0),
               ctrl$verbose,
               inherits(x, "sparseMatrix"),
               package = "quadrupen")
  coefficients <- sparseMatrix(i = out$iA+1,
                               j = out$jA+1,
                               x = c(out$nzeros),
                               dims=c(length(out$lambda1),p))
  active.set   <- sparseMatrix(i = out$iA+1,
                               j = out$jA+1,
                               dims=c(length(out$lambda1),p))
  ## END OF CALL
  if (ctrl$timer) {
    internal.timer <- (proc.time() - cpp.start)[3]
    external.timer <- (proc.time() - r.start)[3]
  } else {
    internal.timer <- NULL
    external.timer <- NULL
  }

  ## ======================================================

  ## BUILDING THE PENFIT OBJECT
  out$converge[out$converge == 0] <- "converged"
  out$converge[out$converge == 1] <- "max # of iterate reached"
  out$converge[out$converge == 2] <- "max # of feature reached"
  monitoring  <- list(it.active      = c(out$it.active ),
                      it.optim       = c(out$it.optim  ),
                      max.grad       = c(out$max.grd   ),
                      status         = c(out$converge  ),
                      pensteps.timer = c(out$timing    ),
                      external.timer = external.timer   ,
                      internal.timer = internal.timer   )

  dimnames(coefficients)[[1]] <- round(c(out$lambda1),3)
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
  residuals <- apply(fitted,2,function(y.hat) y-y.hat)
  r.squared <- 1-colSums(residuals^2)/ifelse(intercept,sum((y-mean(y))^2),sum(y^2))

  return(new("quadrupen",
             coefficients = coefficients   ,
             active.set   = active.set     ,
             intercept    = intercept      ,
             mu           = mu             ,
             normx        = drop(out$normx),
             fitted       = fitted         ,
             residuals    = residuals      ,
             df           = df             ,
             r.squared    = r.squared      ,
             penscale     = penscale       ,
             penalty      = switch(type, "inf" = "l1/linf group Lasso", "l1/l2 group Lasso"),
             naive        = naive          ,
             lambda1      = c(out$lambda1) ,
             lambda2      = lambda2        ,
             monitoring   = monitoring     ,
             control      = ctrl))

}
