DataModel <- R6::R6Class(
  classname = "DataModel",
  private = list(
    names     = NA,
    centered  = NA,
    scaled    = NA
  ),
  public = list(
    ## model-related fields
    X = Matrix(),
    y = numeric(),
    S = Matrix(),
    wx = numeric(),
    wy = numeric(),
    mean_X = numeric(),
    norm_X = double(),
    initialize = 
      function(covariates, outcome, cov_struct,
               intercept=TRUE, standardize=TRUE,
               cov_weights = rep(1,ncol(covariates)),
               obs_weights = rep(1,length(outcome))) {

        ## ===================================================
        ## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
        stopifnot("x has to be of class 'matrix' or 'dgCMatrix'." = 
                  inherits(covariates, c("matrix", "dgCMatrix")))
        stopifnot("NA value in x not allowed." = !any(is.na(covariates)))
        stopifnot("y has to be of type 'numeric'" = is.numeric(outcome))
        stopifnot("x and y have not correct dimensions" = 
                  (nrow(covariates) == length(outcome)))
        stopifnot("struct must be a (square) positive semidefinite matrix." = 
                    all(dim(cov_struct) == p))
        stopifnot("struct must be a (square) positive semidefinite matrix." = 
                    all(eigen(cov_struct, only.values = TRUE)$values >= 0))
        if (!inherits(cov_struct, "sparseMatrix")) 
          cov_struct <- as(cov_struct, "dgCMatrix")
        stopifnot("penscale must have ncol(x) entries" = 
                  (length(cov_weights) == ncol(covariates)))
        stopifnot("covariates weights must be positive" = all(cov_weights > 0))
        stopifnot("observations weights must be positive" = all(obs_weights > 0))
        
        if (is.null(colnames(covariates))) colnames(covariates) <- 1:ncol(covariates)
        
        self$X  <- covariates
        self$y  <- outcome
        self$S  <- cov_struct
        self$wx <- cov_weights
        self$wy <- obs_weights
        private$centered <- intercept
        private$scaled   <- standardize

        ## X and y are normalized but not centered (to keep efficiency with sparse encoding)
        if (intercept) {
          private$names <- c("intercept", colnames(self$X))
          self$mean_X <- colMeans(self$X)
        } else {
          private$names <- colnames(self$X)
          self$mean_X <- rep(0, self$d)
        }
        
        ## normalizing the data
        if (standardize) {
          self$norm_X <-  sqrt(drop(colSums(self$X^2)) - self$n * self$mean_X^2)
          self$X      <- Matrix::colScale(self$X, 1/normx)
          self$mean_X <- self$mean_X/self$norm_X
      } else {
        self$norm_X <- rep(1, self$d)
      }

      self$X      <- Matrix::colScale(self$X, 1/self$wx)
      self$mean_X <- self$mean_X / self$wx

    }
  ), 
  active = list(
    d = function() ncol(self$X),
    n = function() nrow(self$X),
    has_intercept = function() {private$centered},
    is_standardized = function() {private$scaled},
    varnames = function() {private$names}
  )
)

#' @export
GaussianModel <- R6::R6Class(
  classname = "GaussianModel",
  inherit = DataModel,
  public = list(
    mean_y = numeric(),
    norm_y = double(),
    xty    = double(),
    initialize = 
      function(covariates, outcome, cov_struct,
               intercept=TRUE, standardize=TRUE,
               cov_weights = rep(1,ncol(covariates)),
               obs_weights = rep(1,length(outcome))) {
        super$initialize(covariates, outcome, cov_struct, 
                         intercept, standardize, cov_weights, obs_weights)
      if (self$has_intercept) {
        self$mean_y <- mean(self$y)
      } else {
        self$mean_y <- 0
      }
      self$norm_y <- sqrt(sum(self$wy * self$y^2))
      
      if (self$has_intercept) {
        self$xty <- crossprod(self$X, self$y - self$mean_y) - 
          sum(self$y - self$mean_y) * self$mean_X
      } else {
        self$xty <- crossprod(private$x, private$y)
      }
    },
    getL1PenaltyRange = function(length, min_ratio) {
      stopifnot("min.ratio must be non negative." = min_ratio > 0)
      lmax <- max(abs(self$xty))
      lambda <- 10^seq(from=log10(lmax), to=log10(lmax*min_ratio), len=length)  
      lambda
    }
    # loss = function(theta) {
    #   y_hat <- self$X %*% theta
    #   res <- .5 * mean( (self$y - y_hat)^2 )
    #   attr(res, "grad") <- crossprod(self$X, y_hat - self$y)
    #   res
    # }
  ),
  active = list(
    name = function() "Gaussian response (Linear Regression)",
    rss  = function(value) {
      ifelse(private$centered, sum((private$y - mean(private$y))^2), sum(private$y^2))
    }
  )
)

#' 
#' #' @export
#' BinaryModel <- R6::R6Class(
#'   classname = "BinaryModel",
#'   inherit = DataModel,
#'   public = list(
#'     loss = function(theta) {
#'       eta <- self$X %*% theta
#'       res <- sum(log(1 + exp(eta)) - self$y * eta)
#'       attr(res, "grad") <- crossprod(self$X, .sigmoid(eta) - self$y)
#'       res
#'     }
#'   ),
#'   active = list(
#'     name = function() "Binary response (Logistic Regression)"
#'   )
#' )