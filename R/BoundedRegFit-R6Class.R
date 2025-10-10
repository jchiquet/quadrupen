#' @export
#' 
BoundedReg <- R6::R6Class(
  classname = "BoundedReg",
  inherit = QuadrupenFit,
  active  = list(penalty = function(value) "bounded.reg"),
  private = list(naive = NA),
  public  = list(
    initialize =  function(data, naive, lambda1, lambda2) {
      super$initialize(data, lambda1, lambda2)
      private$naive <- naive
      private$data$scaleStruct(lambda2)
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
      private$activeSet <- sparseMatrix(i = out$iB + 1,
                                        j = out$jB + 1,
                                        dims = c(length(private$lambda1),self$ncoef))
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

