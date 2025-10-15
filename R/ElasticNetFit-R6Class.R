#' @export
#' 
ElasticNet <- R6::R6Class(
  classname = "ElasticNet",
  inherit = QuadrupenFit,
  active  = list(penalty = function(value) "elastic.net"),
  private = list(naive = NA),
  public  = list(
    initialize =  function(data, intercept, regParam, naive) {
      super$initialize(data, intercept, regParam)
      private$naive <- naive
      private$optimizer <- elastic_net_cpp
    },
    show = function() {
      super$show()
      if (private$naive) {
        cat("No rescaling of the coefficients (naive Elastic-net).\n")
      } else {
        cat("Coefficients rescaled by (1+lambda2).\n")
      }
    },
    fit = function(control) {

      ## ======================================================
      ## C++ CALL TO ENET_LS
      ## 
      if (control$timer) {cpp.start <- proc.time()}
      out <- private$optimizer(private$data, private$tuning, control)
      timer <- ifelse(control$timer, (proc.time() - cpp.start)[3], NA) 
      ## END OF CALL
      ## ======================================================
      
      private$tuning[[1]] <- out$lambda_l1
      dimnames <- list(round(c(private$tuning[[1]]),3), colnames(private$data$X))
      dims     <- c(length(private$tuning[[1]]),self$ncoef)
      private$mu   <- drop(out$mu)
      private$beta <- sparseMatrix(
          i = out$iA + 1,
          j = out$jA + 1,
          x = c(out$nzeros),
          dims = dims, dimnames = dimnames
        )
      private$activeSet <- sparseMatrix(
        i = out$iA + 1,
        j = out$jA + 1,
        dims = dims, dimnames = dimnames
      )
      private$df <- drop(out$df)
      private$monitoring <- out$monitoring
      private$monitoring$internal.timer <- timer
      private$monitoring$convergence <- 
        sapply(private$monitoring$convergence, status_to_message)
      private$control <- control
    }
  )
)

