#' @export
#' 
ElasticNet <- R6::R6Class(
  classname = "ElasticNet",
  inherit = QuadrupenFit,
  active  = list(
    penalty = function(value) "elastic.net",
    #' @field major_penalty vector of "leading" penalties (either l1 or l2)
    major_penalty = function(value) data.frame(lambda1 = private$lambda1),
    #' @field major_penalty vector of "minor" penalties (either l1 or l2)
    minor_penalty = function(value) setNames(private$lambda2, "lambda2")
  ),
  private = list(
    naive = NA
  ),
  public  = list(
    initialize =  function(data, naive, lambda1, lambda2) {
        super$initialize(data, lambda1, lambda2)
        private$naive <- naive
      },
    show = function() {
      super$show()
      if (private$naive) {
        cat("No rescaling of the coefficients (naive Elastic-net).\n")
      } else {
        cat("Coefficients rescaled by (1+lambda2).\n")
      }
    },
    fit = function(beta0 = NULL, control) {

      if (!is.null(beta0)) {
        stopifnot("beta0 must be a vector with ncol(x) entries." = 
                    (length(beta0) == self$ncoef & is.numeric(beta0)))
      }
      
      D <- Diagonal(x = sqrt(private$lambda2) / sqrt(private$data$wx))
      
      ## ======================================================
      ## STARTING C++ CALL TO ENET_LS
      if (control$timer) {cpp.start <- proc.time()}      
      
      out <- elastic_net2_cpp(
        beta0                     ,
        private$data$X            ,
        private$data$y            ,
        D %*% private$data$S %*% D,
        private$lambda1           ,
        private$data$wx           , 
        private$lambda2           ,
        self$has_intercept        ,
        private$data$xty, private$data$mean_X, private$data$norm_X, private$data$norm_y, 
        private$data$wy      ,
        control$threshold    ,
        control$max.iter     ,
        control$max.feat     ,
        switch(control$method,
               quadra   = 0  ,
               pathwise = 1  ,
               fista    = 2, 0),
        control$verbose      ,
        private$data$sparse_encoding,
        control$usechol,
        control$monitor)
      ## END OF CALL

      if (control$timer) {
        internal.timer <- (proc.time() - cpp.start)[3]
      } else {
        internal.timer <- NULL
      }

      if (!private$naive) {
         out$nzeros <- out$nzeros * (1 + private$lambda2)
         private$mu = private$data$mean_y - (1 + private$lambda2) * drop(out$mu)
       } else {
         private$mu = private$data$mean_y - drop(out$mu);
      }

      private$df      <- drop(out$df)
      private$activeSet <- sparseMatrix(i = out$iA + 1,
                                        j = out$jA + 1,
                                        dims = c(length(private$lambda1),self$ncoef))
      private$beta <- sparseMatrix(i = out$iA + 1,
                                   j = out$jA + 1,
                                   x = c(out$nzeros),
                                   dims = c(length(out$lambda1),self$ncoef),
                                   dimnames = list(round(c(private$lambda1),3),
                                                   colnames(private$data$X)))

      out$converge[out$converge == 0] <- "converged"
      out$converge[out$converge == 1] <- "max # of iterate reached"
      out$converge[out$converge == 2] <- "max # of feature reached"
      out$converge[out$converge == 3] <- "system has become singular"
      private$monitoring <- 
        list(it.active      = c(out$it.active ),
             it.optim       = c(out$it.optim  ),
             max.grad       = c(out$max.grd   ),
             status         = c(out$converge  ),
             pensteps.timer = c(out$timing    ),
             internal.timer = internal.timer   ,
             dist.to.opt    = c(out$delta.hat ),
             dist.to.str    = c(out$delta.star))
    
    }
  )
)
  
