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
    initialize = 
      function(
    x, y, struct, beta0, intercept, normalize, penscale, naive, lambda1, lambda2) {
        super$initialize(
          x, y, struct, beta0, intercept, normalize, penscale, lambda1, lambda2)
        private$naive <- naive
      },
    show = function() {
      super$show()
      if (private$naive) {
        cat("No rescaling of the coefficients (naive).\n")
      } else {
        cat("Coefficients rescaled by (1+lambda2).\n")
      }
    },
    standardize = function() {
      
      if (private$intercept) {
        xbar <- colMeans(private$x)
        ybar <- mean(private$y)
      } else {
        xbar <- rep(0, self$ncoef)
        ybar <- 0;
      }
      
      if (private$normalize) {
        normx <-  sqrt(colSums(private$x^2) - self$nsample * xbar^2)
        private$x <- sweep(private$x, 2, normx, "/")
        xbar <- xbar/normx ;
      } else {
        normx = rep(1, self$ncoef) ;
      }
      normy = sqrt(sum(y^2)) ;
      
      if (any(private$penscale != 1)) {
        private$x <- sweep(private$x, 2, penscale, "/")
        xbar <- xbar/penscale;
      }
      
      if (private$intercept) {
        xty <- crossprod(private$x, private$y - ybar) - sum(private$y - ybar) * xbar
      } else {
        xty <- crossprod(private$x, private$y)
      }
      list(xty = xty, xbar = xbar, normx = normx, normy = normy, ybar = ybar)
      
    },
    fit = function(control) {

      ## ======================================================
      ## STARTING C++ CALL TO ENET_LS
      if (control$timer) {cpp.start <- proc.time()}      
      out <- elastic_net2_cpp(
        private$beta         ,
        private$x            ,
        private$y            ,
        private$struct       ,
        private$lambda1      ,
        private$penscale     ,
        private$lambda2      ,
        private$intercept    ,
        ss$xty, ss$xbar, ss$normx, ss$normy, 
        rep(1,n)             ,
        control$thresh       ,
        control$max.iter     ,
        control$max.feat     ,
        switch(control$method,
               quadra   = 0,
               pathwise = 1,
               fista    = 2, 0),
        control$verbose,
        inherits(x, "sparseMatrix"),
        control$usechol,
        control$monitor)
      ## END OF CALL

      if (control$timer) {
        internal.timer <- (proc.time() - cpp.start)[3]
        external.timer <- (proc.time() - r.start)[3]
      } else {
        internal.timer <- NULL
        external.timer <- NULL
      }
      
      if (!private$naive) {
         out$nonzeros <- out$nonzeros * (1 + private$lambda2)
         private$mu = ss$ybar - (1 + private$lambda2) * drop(out$mu)
       } else {
         private$mu = ss$ybar - drop(out$mu);
      }

      private$df      <- drop(out$df)
      private$normx   <- ss$normx ## drop(out$normx)
      ## private$lambda1 <- out$lambda1
      
      private$activeSet <- sparseMatrix(i = out$iA + 1,
                                        j = out$jA + 1,
                                        dims = c(length(out$lambda1),p))
      private$beta <- sparseMatrix(i = out$iA + 1,
                                   j = out$jA + 1,
                                   x = c(out$nzeros),
                                   dims = c(length(out$lambda1),p))
      
      dimnames(private$beta)[[1]] <- round(c(out$lambda1),3)
      if (is.null(colnames(private$x))) {
        dimnames(private$beta)[[2]] <- 1:p
      } else {
        dimnames(private$beta)[[2]] <- colnames(private$x)
      }

      ## ======================================================
      ## BUILDING THE QUADRUPEN OBJECT
      
      out$converge[out$converge == 0] <- "converged"
      out$converge[out$converge == 1] <- "max # of iterate reached"
      out$converge[out$converge == 2] <- "max # of feature reached"
      out$converge[out$converge == 3] <- "system has become singular"
      private$monitoring  <- list(it.active      = c(out$it.active ),
                          it.optim       = c(out$it.optim  ),
                          max.grad       = c(out$max.grd   ),
                          status         = c(out$converge  ),
                          pensteps.timer = c(out$timing    ),
                          external.timer = external.timer   ,
                          internal.timer = internal.timer   ,
                          dist.to.opt    = c(out$delta.hat ),
                          dist.to.str    = c(out$delta.star))
    
    }
  )
)
  
  