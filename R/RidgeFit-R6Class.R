#' @export
#' 
Ridge <- R6::R6Class(
  classname = "Ridge",
  inherit = QuadrupenFit,
  active  = list(penalty = function(value) "ridge"),
  private = list(
    naive = NA
  ),
  public  = list(
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- ridge_cpp
    },
    fit = function(control) {

      ## ======================================================
      ## C++ CALL TO RIDGE_LS
      ## 
      if (control$timer) {cpp.start <- proc.time()}
      out <- private$optimizer(private$data, private$tuning, control$verbose)
      timer <- ifelse(control$timer, (proc.time() - cpp.start)[3], NA) 
      ## END OF CALL
      ## ======================================================
      
      private$df   <- drop(out$df)
      private$mu   <- drop(out$mu)
      private$beta <- Matrix(
          out$coefficients,
          dimnames = list(round(c(private$tuning[[1]]),3), colnames(private$data$X))
        )
      monitoring  <- list(internal.timer = timer)
    }
  )
)

