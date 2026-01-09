#' Class "ElasticNet"
#' 
#' Class of object returned by the fitting function [elastic.net()]. Inherits fields
#' and methods of [QuadrupenFit]
#' 
#' @seealso [QuadrupenFit], [elastic.net()]
#' 
#' @export
#' 
ElasticNet <- R6::R6Class(
  classname = "ElasticNet",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(penalty = function(value) "elastic.net"),
  public  = list(
    #' @description Initialize a [`ElasticNet`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- elastic_net_cpp
    }
  ),
  private = list(
    rescaled = function() { # Zhou and Hastie Rescaling
      factor <- ( 1 + private$tuning[[2]] )
      list(
        beta = factor * private$beta,
        mu   = factor * private$mu - private$tuning[[2]] * private$data$mean_y
      )
    }
  )
)

#' Class "BoundedRegression"
#' 
#' Class of object returned by the fitting function [bounded.reg()]. Inherits fields
#' and methods of [QuadrupenFit].
#' 
#' @seealso [QuadrupenFit], [bounded.reg()]
#' 
#' @export
BoundedRegression <- R6::R6Class(
  classname = "BoundedRegression",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(penalty = function(value) "bounded.reg"),
  #' @description Initialize a [`BoundedRegression`] model
  #' @param data a [`DataModel`] object
  #' @param intercept a logical; should an intercept be included in the mode?
  #' @param regParam a list with two elements, a vector and a scalar, for the regularization
  public  = list(
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- bounded_reg_cpp
    }
  ),
  private = list(
    rescaled = function() { # Zhou and Hastie Rescaling
      factor <- ( 1 + private$tuning[[2]] )
      list(
        beta = factor * private$beta,
        mu   = factor * private$mu - private$tuning[[2]] * private$data$mean_y
      )
    }
  )
)

#' Class "RidgeRegression"
#' 
#' Class of object returned by the fitting function [ridge()]. Inherits fields
#' and methods of [QuadrupenFit].
#' 
#' @seealso [QuadrupenFit], [ridge()]
#' 
#' @export
RidgeRegression <- R6::R6Class(
  classname = "RidgeRegression",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(penalty = function(value) "ridge"),
  public  = list(
    #' @description Initialize a [`RidgeRegression`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- ridge_cpp
    }
  )
)



#' Class "FusedLasso"
#' 
#' Class of object returned by the fitting function [fused.lasso()]. Inherits fields
#' and methods of [QuadrupenFit].
#' 
#' @seealso [QuadrupenFit], [bounded.reg()]
#' 
#' @export
FusedLasso <- R6::R6Class(
  classname = "BoundedRegression",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(penalty = function(value) "fused.lasso"),
  #' @description Initialize a [`FusedLasso`] model
  #' @param data a [`DataModel`] object
  #' @param intercept a logical; should an intercept be included in the mode?
  #' @param regParam a list with two elements, a vector and a scalar, for the regularization
  public  = list(
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- FusedLasso_cpp
    }
  )
)
