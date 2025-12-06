#' @export
#' 
Ridge <- R6::R6Class(
  classname = "Ridge",
  inherit = QuadrupenFit,
  active  = list(penalty = function(value) "ridge"),
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

      private$tuning[[1]] <- out$tuning_lead
      private$mu          <- drop(out$mu)
      private$beta        <- out$beta
      private$df          <- drop(out$df)
      private$monitoring  <- out$monitoring
      private$monitoring$timer <- timer
      private$control     <- control

      dimnames(private$beta) <- 
        list(round(c(private$tuning[[1]]),3), colnames(private$data$X))
      
    }
  )
)

