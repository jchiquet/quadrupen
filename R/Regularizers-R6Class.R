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
        ifelse(data$sparse_encoding, elastic_net_sparse_cpp, elastic_net_dense_cpp)
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

#' Class "GroupLassoFit"
#' 
#' Class of object returned by the fitting function [group.lasso()]. Inherits fields
#' and methods of [QuadrupenFit]
#' 
#' @seealso [QuadrupenFit], [group.lasso()]
#' 
#' @export
#' 
GroupLassoFit <- R6::R6Class(
  classname = "GroupLassoFit",
  inherit = QuadrupenFit,
  private  = list(group_ = NA, type_ = NA),
  active  = list(
    #' @field penalty character describing the regularizer/penalty
    penalty = function(value) paste0("group.lasso l1/", private$type_),
    #' @field group vector of integers indicating group belonging
    group = function(value) private$group_,
    #' @field type string indicating whether the \eqn{\ell_1/\ell_2}{l1/l2} or the
    #' \eqn{\ell_1/\ell_\infty}{l1/linf} group-Lasso must be fitted.
    type = function(value) private$type_
  ),
  public  = list(
    #' @description Initialize a [`ElasticNetFit`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param group vector of integers indicating group belonging.
    #' @param type string indicating whether the \eqn{\ell_1/\ell_2}{l1/l2} or the
    #' \eqn{\ell_1/\ell_\infty}{l1/linf} group-Lasso must be fitted.
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
    initialize =  function(data, intercept, group, type, regParam) {
      super$initialize(data, intercept, regParam)
      private$group_ <- group
      private$type_  <- type
      private$optimizer <- 
        switch(private$type_,
         "linf" = ifelse(data$sparse_encoding, group_enet_l1linf_sparse_cpp, group_enet_l1linf_dense_cpp),
         "coop" = ifelse(data$sparse_encoding, group_enet_coop_sparse_cpp, group_enet_coop_dense_cpp),
         ifelse(data$sparse_encoding, group_enet_l1l2_sparse_cpp, group_enet_l1l2_dense_cpp)
        )
    },
    #' @description function performing the optimization
    #' @param control list controlling the optimization process
    fit = function(control) {
      ## ======================================================
      ## C++ CALL OPTIMIZER
      ## 
      if (control$timer) {cpp.start <- proc.time()}
      private$stored_fit <- private$optimizer(
        private$data_,
        private$has_intercept_, 
        private$group_,
        private$tuning,
        control
      )
      timer <- ifelse(control$timer, (proc.time() - cpp.start)[3], NA)
      ## END OF CALL
      ## ======================================================
      private$tuning[[1]] <- private$stored_fit$tuning_param[[1]]
      private$intercept_  <- drop(private$stored_fit$intercept)
      private$coef_       <- private$stored_fit$coef
      private$df_         <- drop(private$stored_fit$df)
      private$monitoring  <- private$stored_fit$monitoring
      private$monitoring$timer <- timer
      private$control     <- control
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
  active  = list(
    #' @field penalty character describing the regularizer/penalty
    penalty = function(value) "lava",
    #' @field sparse_coef sparse part of the  decomposition of the coefficients
    sparse_coef = function(value) private$sparse_coef_,
    #' @field dense_coef dense part of the  decomposition of the coefficients
    dense_coef  = function(value) private$coef_- private$sparse_coef_,
    #' @field debias logical, should we rely on the debias coefficient of the regularizer (if available) or not
    debias = function(value) {
      if (missing(value))
        return(private$debias_)
      else {
        stopifnot(is.logical(value))
        private$debias_ <- value
        if (private$debias_) {
          private$coef_ <- private$stored_fit$coef_debias
          private$sparse_coef_ <- private$stored_fit$sparse_coef_debias
          private$intercept_ <- private$stored_fit$intercept_debias
        } else {
          private$coef_ <- private$stored_fit$coef
          private$sparse_coef_ <- private$stored_fit$sparse_coef
          private$intercept_ <- private$stored_fit$intercept
        }
      }
    }
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
    #' @description function performing the optimization
    #' @param control list controlling the optimization process    
    fit = function(control) {
      super$fit(control)
      private$sparse_coef_ <- private$stored_fit$sparse_coef
    }
  )
)
