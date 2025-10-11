#' @export
Ridge <- R6::R6Class(
  classname = "Ridge",
  inherit = QuadrupenFit,
  active  = list(penalty = function(value) "ridge"),
  private = list(
    naive = NA
  ),
  public  = list(
    initialize =  function(lambda) {
      super$initialize(lambda, 0)
      private$data$CholStruct()
    },
    fit = function(data, control) {

      stopifnot("The data object must be an instance of DataModel"
                = inherits(data, "DataModel"))
      private$data <- data
      private$data$scaleStruct(private$lambda)
      
      ## ======================================================
      ## C++ CALL TO RIDGE_LS
      ## 
      if (control$timer) {cpp.start <- proc.time()}
      out <- ridge2_cpp(private$data, private$lambda1, control$verbose)
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

