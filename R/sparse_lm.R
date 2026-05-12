#' Fit a linear model with elastic-net regularization
#'
#' Adjust a linear model with elastic-net regularization, mixing a
#' (possibly weighted) \eqn{\ell_1}{l1}-norm (LASSO) and a
#' (possibly structured) \eqn{\ell_2}{l2}-norm (ridge-like). The
#' solution path is computed at a grid of values for the
#' \eqn{\ell_1}{l1}-penalty, fixing the amount of \eqn{\ell_2}{l2}
#' regularization. See details for the criterion optimized.
#'
#' @inheritParams elastic.net
#' 
#' @param type string indicating the sparse variant to be fitted. 
#' Could be "l1", "mcp" or "scad". Default is "l1".
#' 
#' @param eta real positive scalar for tuning SCAD or MCP penalties. 
#' Default is 3. Ignored when type == "l1".
#'
#' @return an object with class [SparseFit], inheriting from [QuadrupenFit].
#'
#' @seealso See also [QuadrupenFit]
#' 
#' @keywords models, regression
#'
#' @examples
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
#' cor <- 0.75
#' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
#' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
#' Sigma <- Matrix::bdiag(Soo,Sww,Soo,Sww,Soo)
#' diag(Sigma) <- 1
#' n <- 50
#' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
#' y <- 10 + x %*% beta + rnorm(n,0,10)
#'
#' ## Structured Elastic.net without/with an additional l2 regularization term
#' ## and with structuring prior
#' labels <- rep("irrelevant", length(beta))
#' labels[beta != 0] <- "relevant"
#' plot(lasso(x,y), label=labels) ## a mess
#' plot(elastic.net(x,y,lambda2=10), label=labels) ## good guys are selected first
#' plot(elastic.net(x,y,lambda2=10,struct=solve(Sigma)), label=labels) ## even better
#'
#' @export
sparse_lm <- function(x,
                      y,
                      type      = c("l1", "mcp", "scad"),
                      lambda1   = NULL,
                      lambda2   = 0.01,
                      eta       = 3.7,
                      penscale  = rep(1,ncol(x)),
                      struct    = Matrix::Diagonal(ncol(x), 1),
                      intercept = TRUE,
                      normalize = TRUE,
                      refit     = FALSE,
                      nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                      minratio  = ifelse(nrow(x) <= ncol(x), 1e-2, 1e-4),
                      maxfeat   = ifelse(lambda2 < 1e-2, min(nrow(x),ncol(x)), min(4*nrow(x),ncol(x))),
                      beta0     = numeric(ncol(x)),
                      control   = list()) {

  type <- match.arg(type)
  if (type == "mcp") stopifnot(eta > 1)
  if (type == "scad") stopifnot(eta > 2)

  ## ============================================
  ## RECOVER LOW LEVEL CONFIGURATION
  ##
  ctrl <- optim_enet_default(ncol(x))
  ctrl$maxfeat <- maxfeat
  ctrl[names(control)] <- control # default overwritten by user specifications
  ctrl$method <- switch(ctrl$method, quadra = "QUADRA", fista = "FISTA", pgd = "PGD", 0)
  ctrl$factmat <- ctrl$method == "QUADRA"
  ctrl$normalize <- normalize
  ctrl$beta0  <- beta0
  
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
  myModel <- SparseFit$new(
    data      = myData,
    intercept = intercept,
    type      = type,
    regParam  = list(lambda = lambda1, 
                     gamma = lambda2,
                     eta = eta,
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
