#' Class "ElasticNetFit"
#' 
#' Class of object returned by the fitting function [elastic.net()]. Inherits fields
#' and methods of [QuadrupenFit]
#' 
#' @seealso [QuadrupenFit], [elastic.net()]
#' 
#' @export
#' 
ElasticNetFit <- R6::R6Class(
  classname = "ElasticNetFit",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(penalty = function(value) "elastic.net"),
  public  = list(
    #' @description Initialize a [`ElasticNetFit`] model
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
      mean_y <- ifelse(self$has_intercept, mean(private$data$y), 0)
      list(
        beta = factor * private$beta,
        mu   = factor * private$mu - private$tuning[[2]] * mean_y
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
BoundedRegressionFit <- R6::R6Class(
  classname = "BoundedRegressionFit",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(penalty = function(value) "bounded.reg"),
  #' @description Initialize a [`BoundedRegressionFit`] model
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
      mean_y <- ifelse(self$has_intercept, mean(private$data$y), 0)
      list(
        beta = factor * private$beta,
        mu   = factor * private$mu - private$tuning[[2]] * mean_y
      )
    }
  )
)

#' Class "RidgeRegressionFit"
#' 
#' Class of object returned by the fitting function [ridge()]. Inherits fields
#' and methods of [QuadrupenFit].
#' 
#' @seealso [QuadrupenFit], [ridge()]
#' 
#' @export
RidgeRegressionFit <- R6::R6Class(
  classname = "RidgeRegressionFit",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(penalty = function(value) "ridge"),
  public  = list(
    #' @description Initialize a [`RidgeRegressionFit`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- ridge_cpp
    }
  )
)



#' Class "FusedLassoFit"
#' 
#' Class of object returned by the fitting function [fusedlasso()]. Inherits fields
#' and methods of [QuadrupenFit].
#' 
#' @seealso [QuadrupenFit], [bounded.reg()]
#' 
#' @export
FusedLassoFit <- R6::R6Class(
  classname = "FusedLassoFit",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(penalty = function(value) "fused.lasso"),
  #' @description Initialize a [`FusedLassoFit`] model
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
