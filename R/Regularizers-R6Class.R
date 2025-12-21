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
    plot_path = function(xvar = "lambda", log_scale = TRUE,
                         title = paste(self$penalty," path", sep=""),
                         standardize=TRUE, labels = NULL) {
    d <- super$plot_path(xvar, log_scale, title, standardize, labels)
    if (xvar == "lambda") {
      d <- d + xlab(ifelse(log_scale,expression(log[10](lambda[1])),expression(lambda[1])))
      if (log_scale)
          d <- d + scale_x_log10() + annotation_logticks(sides="b")
      } else {
        d <- d + xlab(expression(paste("|",beta[lambda[1]],"|",{}[1]/max[lambda[1]],"|",beta[lambda[1]],"|",{}[1],sep="")))
      }
      d
    }
  )
)

#' @export
#' 
BoundedReg <- R6::R6Class(
  classname = "BoundedReg",
  inherit = QuadrupenFit,
  active  = list(penalty = function(value) "bounded.reg"),
  public  = list(
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- bounded_reg_cpp
    },
    plot_path = function(xvar = "lambda", log_scale = TRUE,
                         title = paste(self$penalty," path", sep=""),
                         standardize=TRUE, labels = NULL) {
      d <- super$plot_path(xvar, log_scale, title, standardize, labels)
      if (xvar == "lambda") {
        d <- d + xlab(ifelse(log_scale,expression(log[10](lambda[infinity])),expression(lambda[infinity])))
        if (log_scale)
          d <- d + scale_x_log10() + annotation_logticks(sides="b")
      } else {
        d <- d + xlab(expression(paste("|",beta[lambda[infinity]],"|",{}[1]/max[lambda[infinity]],"|",beta[lambda[infinity]],"|",{}[1],sep="")))
      }
      d
    }
  )
)

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
    plot_path = function(xvar = "lambda", log_scale = TRUE,
                         title = paste(self$penalty," path", sep=""),
                         standardize=TRUE, labels = NULL) {
      d <- super$plot_path(xvar, log_scale, title, standardize, labels)
      if (xvar == "lambda") {
        d <- d + xlab(ifelse(log_scale,expression(log[10](lambda[2])),expression(lambda[2])))
        if (log_scale)
          d <- d + scale_x_log10() + annotation_logticks(sides="b")
      } else {
        d <- d + xlab(expression(paste("|",beta[lambda[2]],"|",{}[1]/max[lambda[2]],"|",beta[lambda[2]],"|",{}[1],sep="")))
      }
      d
    }
  )
)



