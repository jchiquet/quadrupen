#' Fit a linear model with (sparse) group regularisation (either l1/l2, l1/l-inf or cooperative variant)
#'
#' Adjust a linear model with (sparse) group regularization, that is a
#' mixture of either a (possibly weighted)
#' \eqn{\ell_1/\ell_2}{l1/l2}- or
#' \eqn{\ell_1/\ell_\infty}{l1/linf}-norm, and a (possibly
#' structured) \eqn{\ell_2}{l2}-norm (ridge-like). The solution path
#' is computed at a grid of values for the
#' \eqn{\ell_1/\ell_q}{l1/lq}-penalty. See details for the criterion
#' optimized.
#' 
#' @inheritParams elastic.net
#' 
#' @param alpha real scalar in (0,1); tunes mixture between \eqn{\ell_1}{l1} 
#' group penalties. Default is 0.0 (standard group-lasso).
#'
#' @param group vector of integers indicating group belonging. Must
#' match the number of column in \code{x}. Must be SORTED integers
#' starting from 1.
#'
#' @param type string indicating whether the \eqn{\ell_1/\ell_2}{l1/l2} or the
#' \eqn{\ell_1/\ell_\infty}{l1/linf} group-Lasso must be fitted. Could be "linf" or 
#' "l2", default is "l2"
#'
#' @return an object with class [GroupLassoFit], inheriting from [QuadrupenFit].
#'
#' @seealso See also [QuadrupenFit]
#' 
#' @keywords models, regression
#'
#' @examples
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
#' grp  <- rep(1:5, c(25,10,25,10,25)) 
#' cor <- 0.75
#' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
#' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
#' Sigma <- Matrix::bdiag(Soo,Sww,Soo,Sww,Soo)
#' diag(Sigma) <- 1
#' n <- 50
#' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
#' y <- 10 + x %*% beta + rnorm(n,0,10)
#'
#' ## Group-Lasso without/with an additional l2 regularization term
#' ## and with structuring prior
#' labels <- rep("irrelevant", length(beta))
#' labels[beta != 0] <- "relevant"
#' 
#' \dontrun{
#' ## Standard Group-Lasso
#' plot(group.lasso(x,y,grp), label=labels)
#' plot(group.lasso(x,y,grp, lambda2=.5), label=labels)
#' plot(group.lasso(x,y,grp, lambda2=10), label=labels)
#' plot(group.lasso(x,y,grp, lambda2=10,struct=solve(Sigma)), label=labels)
#' 
#' ## L1/LINF Group-Lasso
#' plot(group.lasso(x, y, grp, type = "linf"), label=labels)
#' plot(group.lasso(x, y, grp, type = "linf", lambda2=.5), label=labels)
#' plot(group.lasso(x, y, grp, type = "linf", lambda2=10), label=labels)
#' plot(group.lasso(x, y, grp, type = "linf", lambda2=10,struct=solve(Sigma)), label=labels)
#' 
#' ## Cooperative-Lasso
#' plot(group.lasso(x, y, grp, type = "coop"), label=labels)
#' plot(group.lasso(x, y, grp, type = "coop", lambda2=.5), label=labels)
#' plot(group.lasso(x, y, grp, type = "coop", lambda2=10), label=labels)
#' plot(group.lasso(x, y, grp, type = "coop", lambda2=10,struct=solve(Sigma)), label=labels)
#' }
#' @export
group.lasso <- function(x,
                        y,
                        group,
                        type      = c("l2", "coop", "linf"),
                        lambda1   = NULL,
                        lambda2   = 0.01,
                        alpha     = 0.0,
                        penscale  = sqrt(tabulate(group)),
                        struct    = Matrix::Diagonal(ncol(x), 1),
                        intercept = TRUE,
                        normalize = TRUE,
                        refit     = FALSE,
                        nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                        minratio  = 1e-2,
                        maxfeat   = ifelse(lambda2 < 1e-2, min(2*nrow(x),ncol(x)), min(4*nrow(x),ncol(x))),
                        beta0     = numeric(ncol(x)),
                        control   = list()) {
  
  ## ============================================
  ## RECOVER LOW LEVEL CONFIGURATION
  ##
  ctrl <- optim_grp_default(ncol(x))
  ctrl$maxfeat <- maxfeat
  if (!is.null(control$method)) if (control$method != "quadra") ctrl$threshold <- 1e-2
  ctrl[names(control)] <- control # default overwritten by user specifications
  ctrl$method  <- switch(ctrl$method, quadra = "QUADRA", fista = "FISTA", pgd = "PGD", 0)
  ctrl$usechol <- FALSE
  ctrl$normalize <- normalize
  ctrl$beta0  <- beta0
  stopifnot(alpha < 1 && alpha >= 0)
  stopifnot(!is.unsorted(group))
  
  ## ============================================
  ## INSTANTIATE THE DATA MODEL
  ##
  myData <- DataModel$new(
    covariates  = x,
    outcome     = y,
    cov_struct  = struct
  )
  
  ## ============================================
  ## INSTANTIATE THE PENALIZED MODEL
  ##
  myModel <- GroupLassoFit$new(
    data      = myData,
    intercept = intercept,
    group     = group,
    type      = match.arg(type),
    regParam  = list(lambda = lambda1, 
                     gamma  = lambda2,
                     alpha  = alpha,
                     lambda_factor = penscale, 
                     min_ratio = minratio, n_lambda = nlambda1)
  )
  
  ## ============================================
  ## FIT THE MODEL WITH ACTIVE SET ALGORITHM
  ##
  if (ctrl$verbose > 0) cat("\nModel fitting and optimization")
  myModel$fit(ctrl)
  
  ## ============================================
  ## POSTREATMENT + SEND BACK THE RESULTING MODEL
  ##
  if (ctrl$verbose > 0) cat("\nPost-treatment")
  myModel$debias <- refit
  myModel$criteria()
  myModel
}
