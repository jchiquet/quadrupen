#' @export
#' 
ElasticNet <- R6::R6Class(
  classname = "ElasticNet",
  inherit = QuadrupenFit,
  active  = list(penalty = function(value) "elastic.net"),
  public  = list(
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- elastic_net_cpp
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

      private$tuning[[1]] <- out$tuning_lead
      private$mu          <- drop(out$mu)
      private$beta        <- out$beta
      private$df          <- drop(out$df)
      private$monitoring  <- out$monitoring
      private$monitoring$timer <- timer
      private$control     <- control

      private$monitoring$convergence <- 
        sapply(private$monitoring$convergence, status_to_message)
      dimnames(private$beta) <- 
        list(round(c(private$tuning[[1]]),3), colnames(private$data$X))
      
    }
  )
)

