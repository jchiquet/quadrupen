#' @export
#' 
BoundedReg <- R6::R6Class(
  classname = "BoundedReg",
  inherit = QuadrupenFit,
  active  = list(penalty = function(value) "bounded.reg"),
  private = list(naive = NA),
  public  = list(
    initialize =  function(data, intercept, lambda1, lambda2, naive) {
      super$initialize(data, intercept, lambda1, lambda2)
      private$naive <- naive
    },
    show = function() {
      super$show()
      if (private$naive) {
        cat("No rescaling of the coefficients (naive bounded regression).\n")
      } else {
        cat("Coefficients rescaled by (1+lambda2).\n")
      }
    },
    fit = function(control) {

      ## ======================================================
      ## C++ CALL TO BREG_LS
      ## 
      if (control$timer) {cpp.start <- proc.time()}
      out <- 
        bounded_reg_cpp(
          private$data, 
          private$lambda1,
          private$lambda2,
          control
        )
      timer <- ifelse(control$timer, (proc.time() - cpp.start)[3], NA) 
      ## END OF CALL
      ## ======================================================

      private$df      <- drop(out$df)
      private$lambda1 <- out$lambda1
      private$activeSet <- sparseMatrix(i = out$iB + 1,
                                        j = out$jB + 1,
                                        dims = c(length(out$lambda1),self$ncoef))
      private$mu   <- drop(out$mu)
      private$beta <- Matrix(out$coef, 
            dimnames = list(round(c(private$lambda1),3), colnames(private$data$X)))
      
      private$monitoring <- out$monitoring
      private$monitoring$internal.timer <- timer
      private$monitoring$convergence <- 
        sapply(private$monitoring$convergence, status_to_message)
      private$control <- control
    }
  )
)

