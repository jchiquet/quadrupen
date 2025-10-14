#' @export
Ridge <- R6::R6Class(
  classname = "Ridge",
  inherit = QuadrupenFit,
  active  = list(penalty = function(value) "ridge"),
  private = list(
    naive = NA
  ),
  public  = list(
    initialize =  function(data, intercept, lambda) {
      super$initialize(data, intercept, lambda, 0)
    },
    fit = function(control) {

      ## ======================================================
      ## C++ CALL TO RIDGE_LS
      ## 
      if (control$timer) {cpp.start <- proc.time()}
      out <- ridge_cpp(private$data, private$lambda1, control$verbose)
      timer <- ifelse(control$timer, (proc.time() - cpp.start)[3], NA) 
      ## END OF CALL
      ## ======================================================
      
      private$df   <- drop(out$df)
      private$mu   <- drop(out$mu)
      private$beta <- Matrix(out$coefficients,
                             dimnames = list(round(c(private$lambda1),3),
                                        colnames(private$data$X)))
      monitoring  <- list(internal.timer = timer)
    }
  )
)

