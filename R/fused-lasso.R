#' A function for fitting generalized fused-Lasso problems
#' 
#' This function fits the standard version of the fused lasso. 
#' It can take a general matrix x and provides for possible weights on the 
#' `lambda1` and `lambda2` penalties. 
#' 
#' @inheritParams elastic.net
#' 
#' @param struct Description of the graph that corresponds to the `lambda2` penalty structure. 
#' If \code{NULL} (the default) a chain graph is assumed, like in the standard fused-lasso.
#' If a matrix is given, interpreted as 
#' a symmetric adjacency matrix
#' @param pen_fused penalty used for fusing the variables (either L1, L2 or Huber). Default is L1 
#' 
#' @param control list of argument controlling low level options of
#' the algorithm:
#' * `verbose`: logical; ener verbose mode
#' * `timer`: logical; use to record the timing of the
#' algorithm. Default is `FALSE`.
#' * `maxiterin` Maximum number of iterations in the inner loop to run.
#' * `maxiterout` Maximum number of iterations in the outer loop to run.
#' * `maxactivation` Maximum number of previously inactive variables to activate at the same time
#' * `accuracy` Accuracy at which the algorithm will stop.
#' @param fusionCheck Should the fused sets be checked for breaking up? 
#' @param verbose Should the function give some output what it is doing?
#' 
#' @returns A list with the elements
#' * beta A sparse matrix of type \code{dgCMatrix} that returns the solutions for each value 
#' of \code{lambda1} and \code{lambda2}. Each column is a separate solution.
#' * lambda1 The vector of the lambda1 values that were used.
#' * lambda2 The vector of the lambda2 values that were used.
#' * success A logical vector inidicating if the algorithm converged.
#' * Intercept The intercept of the model; 0 if no intercept was included
#' 
#' @author Holger Hoefling, refactoring by Julien C.
#' @seealso \code{\link{testData}}
#' 
#' @examples 
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
#' res1 <- fusedlasso(x, y, lambda2=0.1)
#' G1 = quadrupen:::chain_graph(ncol(x))
#' res2 <- fusedlasso(x, y, lambda2=1e-1, graph = G1)
#' G2 <- igraph::make_ring(ncol(x)) |> igraph::as_adjacency_matrix(sparse = FALSE)
#' res3 <- fusedlasso(x, y, lambda2=0.1, graph = G2)
#' matplot(log10(res1$tuning_param$l1), t(as.matrix(res1$beta)), type="l")
#' matplot(log10(res2$tuning_param$l1), t(as.matrix(res2$beta)), type="l")
#' matplot(log10(res3$tuning_param$l1), t(as.matrix(res3$beta)), type="l")
#' 
#' @keywords models regression multivariate
#' 
#' @importFrom methods as
#' @export 
fusedlasso <- function(x, 
                       y, 
                       lambda1   = NULL,                       
                       lambda2   = 0.01,
                       pen_fused  =  c("L1", "L2", "Huber"),
                       penscale  = rep(1,ncol(x)), ## doesnt work at the moment
                       struct    = NULL,
                       intercept = TRUE,
                       normalize = TRUE,
                       debiasing = c("none", "standard", "ridge"),
                       nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                       minratio  = ifelse(nrow(x) <= ncol(x), 1e-3, 1e-4),
                       maxfeat   = ifelse(lambda2 < 1e-2, min(nrow(x),ncol(x)), min(2*nrow(x),ncol(x))),
                       beta0     = rep(0, ncol(x)),
                       control   = list()) {

  ## ============================================
  ## RECOVER LOW LEVEL CONFIGURATION
  ##
  ctrl <- ctrl_fused_default(ncol(x))
  ctrl$maxfeat <- maxfeat
  ctrl[names(control)] <- control # default overwritten by user specifications
  ctrl$beta0  <- beta0
  ctrl$mu0    <- 0.001
  ctrl$intercept <- intercept
  ctrl$pen_fused  <- match.arg(pen_fused)
  ctrl$rescaling <- match.arg(debiasing)
  ctrl$fusioncheck <- 
    switch(ctrl$fusioncheck, "all" = 2L, "active" = 1L, "none" = 0L, "naive" = -1L, -2L)
  if (ctrl$pen_fused == "L2")  ctrl$fusioncheck <- 0
  if (is.null(struct))   struct <- chain_graph(ncol(x))
  if (is.matrix(struct)) struct <- conn_from_adj(struct)
  
  ## ============================================
  ## INSTANTIATE THE DATA MODEL
  ##
  if (ctrl$verbose > 0) cat ("\nData pretreatment")
  myData <- GaussianModel$new(
    covariates  = as(x, "dgCMatrix"),
    outcome     = y,
    cov_struct  = struct,
    cov_weights = penscale
  )
  myData$standardize(intercept, normalize)
  myData$getSufficientStat()

  ## ===========================
  ## Tuning parameters
  if (is.null(lambda1)) {
    ctrl$adjust <- TRUE
    lambda1 <- c(1, exp((1:(nlambda1-1)) * log(minratio)/(nlambda1-1)))
  }
  if (length(lambda2) != 1 & length(lambda1) != length(lambda2)) {
    stop("Give either 1 value for lambda2 or same number as lambda2")
  }
  if (any(lambda2 <= 0)) stop("Lambda2 has to have only positive elements")
  if (length(lambda2) == 1) lambda2 <- rep(lambda2, length(lambda1))
  tuningParam <- list("l1" = nrow(x) * lambda1, "l2" = nrow(x) * lambda2)  

  ## ===========================
  ## Call to the main function
  ## 
  res <- FusedLasso_cpp(myData, tuningParam, ctrl)
  
  ## ===========================
  ## Post-Treatments
  ## 
  if (intercept) {
    res$mu   <- res$beta[nrow(res$beta)]
    res$beta <- res$beta[1:(nrow(res$beta) - 1) , , drop = FALSE]
  } else {
    res$mu <- rep(0, length(lambda1))
  }
  res$beta <- res$beta / myData$norm_X
  res$tuning_param <- lapply(res$tuning_param, function(param) param/nrow(x))

  res
}
