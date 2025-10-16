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
    C = matrix(),
    S = Matrix(),
    wx = numeric(),
    wy = numeric(),
    mean_X = numeric(),
    norm_X = double(),
    mean_y = numeric(),
    norm_y = double(),
    initialize = 
      function(covariates, outcome, cov_struct,
               cov_weights = rep(1,ncol(covariates)),
               obs_weights = rep(1,length(outcome)), check_args = TRUE) {
        
        ## ===================================================
        ## CHECKS TO (PARTIALLY) AVOID CRASHES OF THE C++ CODE
        ##
        if (check_args) {
          stopifnot("x has to be of class 'matrix' or 'dgCMatrix'." = 
                      inherits(covariates, c("matrix", "dgCMatrix")))
          stopifnot("NA value in x not allowed." = !any(is.na(covariates)))
          stopifnot("y has to be of type 'numeric'" = is.numeric(outcome))
          stopifnot("x and y have not correct dimensions" = 
                      (nrow(covariates) == length(outcome)))
          if (!is.null(cov_struct)) {
            stopifnot("struct must be a (square) positive semidefinite matrix." = 
                        all(dim(cov_struct) == ncol(covariates)))
            stopifnot("struct must be a (square) positive semidefinite matrix." = 
                        all(eigen(cov_struct, only.values = TRUE)$values >= 0))
            if (!inherits(cov_struct, "sparseMatrix")) 
              cov_struct <- as(cov_struct, "CsparseMatrix")
          }
          stopifnot("penscale must have ncol(x) entries" = 
                      (length(cov_weights) == ncol(covariates)))
          stopifnot("covariates weights must be positive" = all(cov_weights > 0))
          stopifnot("observations weights must be positive" = all(obs_weights > 0))
        } 
        ## ===================================================
        
        if (is.null(colnames(covariates))) colnames(covariates) <- 1:ncol(covariates)
        self$X  <- covariates
        self$y  <- outcome
        self$S  <- cov_struct
        self$wx <- cov_weights
        self$wy <- obs_weights
        private$centered <- FALSE
        private$scaled   <- FALSE
      },
    standardize = function(intercept, normalize) {
      ## X and y are not centered to keep efficiency with sparse encoding

      private$centered <- intercept
      private$scaled   <- normalize
      
      if (intercept) {
        private$names <- c("intercept", colnames(self$X))
        self$mean_X <- colMeans(self$X)
        self$mean_y <- mean(self$y)
      } else {
        private$names <- colnames(self$X)
        self$mean_X <- rep(0, self$d)
        self$mean_y <- 0
      }
      
      ## normalizing the data
      self$norm_y <- sqrt(sum(self$wy * self$y^2))
      if (normalize) {
        self$norm_X <-  sqrt(drop(colSums(self$X^2)) - self$n * self$mean_X^2)
      } else {
        self$norm_X <- rep(1, self$d)
      }

      self$X      <- Matrix::colScale(self$X, 1/(self$norm_X * self$wx))
      self$mean_X <- self$mean_X / (self$wx * self$norm_X)
      if (!inherits(self$X, "sparseMatrix")) self$X <- as.matrix(self$X)
      ##
      ## ===================================================
    },
    scaleStruct = function(lambda) {
      self$S <- dimScale(self$S, sqrt(lambda/self$wx))
    },
    CholStruct = function() {
      self$C <- as.matrix(chol(self$S))
    }
  ), 
  active = list(
    d = function() ncol(self$X),
    n = function() nrow(self$X),
    has_intercept = function() {private$centered},
    is_centered = function() {private$centered},
    is_scaled = function() {private$scaled},
    sparse_encoding = function() {inherits(private$X, "sparseMatrix")},
    varnames = function() {private$names}
  )
)

#' @export
GaussianModel <- R6::R6Class(
  classname = "GaussianModel",
  inherit = DataModel,
  public = list(
    xty    = double(),
    getSufficientStat = function() {
      self$xty <- as.numeric(crossprod(self$X, self$y - self$mean_y) - 
                               sum(self$y - self$mean_y) * self$mean_X)
    },
    show = function() {
      cat("Gaussian Data," ,
          ifelse(private$centered, "centered", "not centered"), "and",
          ifelse(private$scaled, "scaled", "not scaled.\n"))
      invisible(self)
    },
    #' @description User friendly print method
    print = function() { self$show() },
    splitTrainTest = function(
    nfolds = 10,
    folds  = split(sample(1:self$n), rep(1:nfolds, length = self$n))
    ) {
      ## un-normalize data 
      X <- Matrix::colScale(self$X, self$norm_X * self$wx)
      if (!inherits(X, "sparseMatrix")) X <- as.matrix(X)
      ## create the list of split each compose with couple of Train/Test
      lapply(folds, function(omit) {
        trainData <- GaussianModel$new(X[-omit, ], self$y[-omit], self$S, self$wx, self$wy[-omit])
        testData  <- GaussianModel$new(X[omit, ], self$y[omit], self$S, self$wx, self$wy[omit])
        list(trainData = trainData, testData = testData,
             trainID = setdiff(1:self$n, omit), testID = omit)
      })
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
    rss  = function(value) {sum((self$y - self$mean_y)^2)}
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