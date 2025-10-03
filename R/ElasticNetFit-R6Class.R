#' @export
#' 
ElasticNet <- R6::R6Class(
  classname = "ElasticNet",
  inherit = QuadrupenFit,
  active  = list(
    penalty = function(value) "elastic.net",
    #' @field major_penalty vector of "leading" penalties (either l1 or l2)
    major_penalty = function(value) data.frame(lambda1 = private$lambda1),
    #' @field major_penalty vector of "minor" penalties (either l1 or l2)
    minor_penalty = function(value) setNames(private$lambda2, "lambda2")
  ),
  private = list(
    naive = NA
  ),
  public  = list(
    initialize =  function(data, naive, lambda1, lambda2) {
        super$initialize(data, lambda1, lambda2)
        private$naive <- naive
      },
    show = function() {
      super$show()
      self$data$scaleStruct(lambda2)
      if (private$naive) {
        cat("No rescaling of the coefficients (naive Elastic-net).\n")
      } else {
        cat("Coefficients rescaled by (1+lambda2).\n")
      }
    },
    fit = function(beta0 = NULL, control) {

      if (!is.null(beta0)) {
        stopifnot("beta0 must be a vector with ncol(x) entries." = 
                    (length(beta0) == self$ncoef & is.numeric(beta0)))
      }

      ## ======================================================
      ## C++ CALL TO ENET_LS
      ## 
      if (control$timer) {cpp.start <- proc.time()}
      out <- 
        elastic_net2_cpp(
          beta0, 
          private$data, 
          private$lambda1,
          private$lambda2,
          control
        )
      timer <- ifelse(control$timer, (proc.time() - cpp.start)[3], NA) 
      ## END OF CALL
      ## ======================================================

      private$df      <- drop(out$df)
      private$activeSet <- sparseMatrix(i = out$iA + 1,
                                        j = out$jA + 1,
                                        dims = c(length(private$lambda1),self$ncoef))
      private$mu   <- drop(out$mu)
      private$beta <- sparseMatrix(i = out$iA + 1,
                                   j = out$jA + 1,
                                   x = c(out$nzeros),
                                   dims = c(length(out$lambda1),self$ncoef),
                                   dimnames = list(round(c(private$lambda1),3),
                                                   colnames(private$data$X)))
      private$monitoring <- out$monitoring
      private$monitoring$internal.timer <- timer
      private$monitoring$convergence <- 
        sapply(private$monitoring$convergence, status_to_message)

    }
  )
)
  
