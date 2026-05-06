#' @rdname group_sparse_lm
#' @export
group_lasso <- 
  function(x,
           y,
           group,
           lambda1   = NULL,
           penscale  = sqrt(tabulate(group)),
           intercept = TRUE,
           normalize = TRUE,
           nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
           minratio  = 1e-2,
           maxfeat   = ncol(x),
           beta0     = numeric(ncol(x)),
           control   = list(method = "quadra")) {
    
    out <- group_sparse_lm(
      x,
      y,
      group,
      type      = "l2",
      lambda1   = lambda1,
      lambda2   = 0.0,
      alpha     = 0.0,
      penscale  = penscale,
      struct    = Matrix::Diagonal(ncol(x), 1),
      intercept = intercept,
      normalize = normalize,
      refit     = FALSE,
      nlambda1  = nlambda1,
      minratio  = minratio,
      maxfeat   = maxfeat,
      beta0     = beta0,
      control   = control)
    
    out    
  }

#' @rdname group_sparse_lm
#' @export
group_l1linf <- 
  function(x,
           y,
           group,
           lambda1   = NULL,
           penscale  = sqrt(tabulate(group)),
           intercept = TRUE,
           normalize = TRUE,
           nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
           minratio  = 1e-2,
           maxfeat   = ncol(x),
           beta0     = numeric(ncol(x)),
           control   = list()) {
    
    out <- group_sparse_lm(
      x,
      y,
      group,
      type      = "linf",
      lambda1   = lambda1,
      lambda2   = 0.0,
      alpha     = 0.0,
      penscale  = penscale,
      struct    = Matrix::Diagonal(ncol(x), 1),
      intercept = intercept,
      normalize = normalize,
      refit     = FALSE,
      nlambda1  = nlambda1,
      minratio  = minratio,
      maxfeat   = maxfeat,
      beta0     = beta0,
      control   = control)
    
    out    
  }

#' @rdname group_sparse_lm
#' @export
coop_lasso <- 
  function(x,
           y,
           group,
           lambda1   = NULL,
           penscale  = sqrt(tabulate(group)),
           intercept = TRUE,
           normalize = TRUE,
           nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
           minratio  = 1e-2,
           maxfeat   = ncol(x),
           beta0     = numeric(ncol(x)),
           control   = list()) {
    
    out <- group_sparse_lm(
      x,
      y,
      group,
      type      = "coop",
      lambda1   = lambda1,
      lambda2   = 0.0,
      alpha     = 0.0,
      penscale  = penscale,
      struct    = Matrix::Diagonal(ncol(x), 1),
      intercept = intercept,
      normalize = normalize,
      refit     = FALSE,
      nlambda1  = nlambda1,
      minratio  = minratio,
      maxfeat   = maxfeat,
      beta0     = beta0,
      control   = control)
    
    out    
  }

#' @rdname group_sparse_lm
#' @export
sparse_group_lasso <- 
  function(x,
           y,
           group,
           lambda1   = NULL,
           alpha     = 0.5,
           penscale  = sqrt(tabulate(group)),
           intercept = TRUE,
           normalize = TRUE,
           nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
           minratio  = 1e-2,
           maxfeat   = ncol(x),
           beta0     = numeric(ncol(x)),
           control   = list()) {
    
    out <- group_sparse_lm(
      x,
      y,
      group,
      type      = "l2",
      lambda1   = lambda1,
      lambda2   = 0.0,
      alpha     = alpha,
      penscale  = penscale,
      struct    = Matrix::Diagonal(ncol(x), 1),
      intercept = intercept,
      normalize = normalize,
      refit     = FALSE,
      nlambda1  = nlambda1,
      minratio  = minratio,
      maxfeat   = maxfeat,
      beta0     = beta0,
      control   = control)
    
    out    
  }

#' @rdname group_sparse_lm
#' @export
sparse_group_l1linf <- 
  function(x,
           y,
           group,
           lambda1   = NULL,
           alpha     = 0.0,
           penscale  = sqrt(tabulate(group)),
           intercept = TRUE,
           normalize = TRUE,
           nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
           minratio  = 1e-2,
           maxfeat   = ncol(x),
           beta0     = numeric(ncol(x)),
           control   = list()) {
    
    out <- group_sparse_lm(
      x,
      y,
      group,
      type      = "linf",
      lambda1   = lambda1,
      lambda2   = 0.0,
      alpha     = alpha,
      penscale  = penscale,
      struct    = Matrix::Diagonal(ncol(x), 1),
      intercept = intercept,
      normalize = normalize,
      refit     = FALSE,
      nlambda1  = nlambda1,
      minratio  = minratio,
      maxfeat   = maxfeat,
      beta0     = beta0,
      control   = control)
    
    out    
  }

#' @rdname group_sparse_lm
#' @export
sparse_coop_lasso <- 
  function(x,
           y,
           group,
           lambda1   = NULL,
           alpha     = 0.0,
           penscale  = sqrt(tabulate(group)),
           intercept = TRUE,
           normalize = TRUE,
           nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
           minratio  = 1e-2,
           maxfeat   = ncol(x),
           beta0     = numeric(ncol(x)),
           control   = list()) {
    
    out <- group_sparse_lm(
      x,
      y,
      group,
      type      = "coop",
      lambda1   = lambda1,
      lambda2   = 0.0,
      alpha     = alpha,
      penscale  = penscale,
      struct    = Matrix::Diagonal(ncol(x), 1),
      intercept = intercept,
      normalize = normalize,
      refit     = FALSE,
      nlambda1  = nlambda1,
      minratio  = minratio,
      maxfeat   = maxfeat,
      beta0     = beta0,
      control   = control)
    
    out    
  }
