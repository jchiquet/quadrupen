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
      private$optimizer <- 
        ifelse(data$sparse_encoding,
               elastic_net_sparse_cpp,
               elastic_net_dense_cpp)
    }
  ),
  private = list(
    rescaled = function(beta = private$beta, mu = private$mu) { # Zhou and Hastie Rescaling
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
      private$optimizer <- bounded_regression_cpp
    },
    #' @description function performing the optimization
    #' @param control list controlling the optimization process
    fit = function(control) {
      super$fit(control)
      private$beta <-  do.call(rbind, private$beta)
    }
  ),
  private = list(
    rescaled = function(beta = private$beta, mu = private$mu) { # Zhou and Hastie Rescaling
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
    },
    #' @description function performing the optimization
    #' @param control list controlling the optimization process
    fit = function(control) {
      super$fit(control)
      private$beta <-  do.call(rbind, private$beta)
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

#' Class "LavaFit"
#' 
#' Class of object returned by the fitting function [lava()]. Inherits fields
#' and methods of [QuadrupenFit]
#' 
#' @seealso [QuadrupenFit], [lava()]
#' 
#' @export
#' 
LavaFit <- R6::R6Class(
  classname = "LavaFit",
  inherit = QuadrupenFit,
  #' @field penalty character describing the regularizer/penalty
  active  = list(
    penalty = function(value) "lava",
    sparse_coef = function(value) sparse_coef_,
    dense_coef  = function(value) dense_coef_
    ),
  private = list(sparse_coef_ = NA, dense_coef_ = NA),
  public  = list(
    #' @description Initialize a [`LavaFit`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, regParam) {
      super$initialize(data, intercept, regParam)
      private$optimizer <- lava_dense_cpp
    },
    fit = function(control) {
      out <- super$fit(control)
      private$sparse_coef_ <- out$delta
      private$dense_coef_  <-  do.call(rbind, private$beta)
      private$beta <- private$sparse_coef_ + private$dense_coef_
    }
  )
)
