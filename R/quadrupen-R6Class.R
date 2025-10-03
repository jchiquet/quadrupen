#' Class "QuadrupenFit"
#'
#' Class of object returned by any fitting function of the
#' \pkg{quadrupen} package (\code{elastic.net} or
#' \code{bounded.reg}).
#' 
#' This class comes with the usual \code{predict(object, newx, ...)},
#' \code{fitted(object, ...)}, \code{residuals(object, ...)},
#' \code{print(object, ...)}, \code{show(object)} and
#' \code{deviance(object, ...)} generic (undocumented) methods.
#'
#' A specific plotting method is available and documented
#' (\code{\link{plot.quadrupen}}).
#'
#' @param coefficients Matrix (class \code{"dgCMatrix"}) of
#' coefficients with respect to the original input. The number of
#' rows corresponds the length of \code{lambda1}.
#'
#' @param activeSet Matrix (class \code{"dgCMatrix"}, generally
#' sparse) indicating the 'active' variables, in the sense that they
#' activate the constraints. For the \code{\link{elastic.net}}, it
#' corresponds to the nonzero variables; for the
#' \code{\link{bounded.reg}} function, it is the set of variables
#' reaching the boundary along the path of solutions.
#'
#' @param intercept logical; indicates if an intercept has
#'  been included to the model.
#'
#' @param mu A vector (class \code{"numeric"})
#' containing the successive values of the (unpenalized) intercept.
#' Equals to zero if \code{intercept} has been set to \code{FALSE}.
#'
#' @param normx Vector (class \code{"numeric"}) containing the
#' square root of the sum of squares of each column of the design
#' matrix.
#'
#' @param penscale Vector \code{"numeric"} with real positive
#' values that have been used to weight the penalty tuned by
#' \eqn{\lambda_1}{lambda1}.
#'
#' @param penalty Object of class \code{"character"}
#' indicating the method used (\code{"elastic-net"} or \code{"bounded
#' regression"}).
#'
#' @param naive logical; was the \code{naive} mode on?
#'
#' @param lambda1 Vector (class \code{"numeric"}) of penalty
#' levels (either \eqn{\ell_1}{l1} or \eqn{\ell_\infty}{l-infinity})
#' for which the model has eventually been fitted.
#'
#' @param lambda2 Scalar (class \code{"numeric"}) for the
#' amount of \eqn{\ell_2}{l2} (ridge-like) penalty.
#'
#' @param control Object of class \code{"list"} with low
#' level options used for optimization.
#'
#' @param monitoring List (class \code{"list"}) which
#' contains various indicators dealing with the optimization
#' process.
#'
#' @param residuals Matrix of residuals, each column
#' corresponding to a value of \code{lambda1}.
#'
#' @param df Estimated degree of freedoms for the successive
#' \code{lambda1}.  Only available for 'elastic.net' using tCholesky
#' factorization.
#'
#' @param r.squared Vector (class \code{"numeric"}) given the
#' coefficient of determination as a function of lambda1.
#'
#' @param fitted Matrix of fitted values, each column
#' corresponding to a value of \code{lambda1}.  
#'
#' @seealso See also \code{\link{plot.quadrupen}}.
#'
#' @importFrom stats fitted predict residuals deviance
#' @export
#' 
QuadrupenFit <- R6Class(
  classname = "QuadrupenFit",
  ## ____________________________________________________
  ## 
  ## PRIVATE MEMBERS
  ## ____________________________________________________
  private = list(
    data        = NA,
    beta        = Matrix()  ,
    mu          = numeric() ,
    activeSet   = Matrix()  ,
    df          = numeric() ,
    lambda1     = numeric() ,
    lambda2     = numeric() ,
    control     = list()    ,
    monitoring  = list()
  ),
  ## ____________________________________________________
  ## 
  ## ACTIVE BINDINGS MEMBERS
  ## ____________________________________________________
  active = list(
    ncoef = function(value) {private$data$d},
    nsample = function(value) {private$data$n},
    has_intercept = function(value) {private$data$has_intercept},
    is_standardized = function(value) {private$data$is_standardized},
    fitted = function(value) {
      if (self$has_intercept) {
### TODO - normalize x back        
        res <- sweep(tcrossprod(private$data$X, private$beta),2L,-private$mu,check.margin=FALSE)
      } else {
        private$mu <- 0
        res <- tcrossprod(private$data$X, private$beta)
      }
      res
    },
    coefficients = function(value) {private$beta},
    residuals = function(value) {apply(self$fitted, 2, function(y_hat) private$data$y - y_hat)},
    deviance = function(value) {colSums(self$residuals^2)},
    degrees_freedom = function(value) {
      private$df + ifelse(private$intercept, 1L, 0L)
    },
### TODO - only valid for Gaussian models
    r_squared = function(value) {
      1 - colSums(self$residuals^2) / private$data$rss
    }
  ),

  ## ____________________________________________________
  ## 
  ## PUBLIC MEMBERS
  ## ____________________________________________________
  public  = list(
    initialize = function(data, lambda1, lambda2) {

      ## ===================================================
      ## CHECKS TO AVOID CRASHES OF THE C++ CODE
      stopifnot("The data object must be an instance of DataModel"
                = inherits(data, "DataModel"))
      stopifnot("entries inlambda1 must all be postive." =
                  (all(lambda1 > 0)))
      stopifnot("lambda1 values must be sorted in decreasing order." =
                  !is.unsorted(rev(lambda1)))
      stopifnot("lambda2 must be a scalar." = 
                  (length(lambda2) == 1 & inherits(lambda2, "numeric")))
      stopifnot("lambda2 must be a non negative scalar." =  
                  (lambda2 >= 0))
      private$data     <- data
      private$lambda1  <- lambda1
      private$lambda2  <- lambda2
    },
    show = function() {
      cat("Linear regression with", self$penalty, "penalizer.\n")
      if (self$has_intercept) {
        cat("- number of coefficients:", self$ncoef,"+ intercept\n")
      } else {
        cat("- number of coefficients:", self$ncoef,"(no intercept)\n")
      }
      cat("- penalty parameter ",names(self$major_penalty), ": ",
          length(self$major_penalty), " points from ",
          format(max(self$major_penalty), digits = 3)," to ",
          format(min(self$major_penalty), digits = 3),"\n", sep="")
      cat("- penalty parameter ",names(self$minor_penalty),": ", self$minor_penalty, "\n", sep="")
      invisible(self)
    },
    #' @description User friendly print method
    print = function() { self$show() },
    predict = function(newx = NULL, ... ) {
      if (is.null(newx)) {
        res <- self$fitted
      } else {
        res <- sweep(newx %*% t(private$beta),2L,-private$mu, check.margin=FALSE)
      }
      res
    },
    #' Plot method for a quadrupen object
    #'
    #' Produce a plot of the solution path of a \code{quadrupen} fit.
    #'
    #' @usage plot.quadrupen(x, y, xvar = "lambda",
    #'         main = self$penalty," path", sep=""),
    #'         log.scale = TRUE, standardize=TRUE, reverse=FALSE,
    #'         labels = NULL, plot = TRUE, ...)
    #' @param x output of a fitting procedure of the \pkg{quadrupen}
    #' package (\code{\link{elastic.net}} or \code{\link{bounded.reg}}
    #' for the moment). Must be of class \code{quadrupen}.
    #' @param y used for S4 compatibility.
    #' @param xvar variable to plot on the X-axis: either \code{"lambda"}
    #' (\eqn{\lambda_1}{lambda1} penalty level or
    #' \eqn{\lambda_2}{lambda2} for ridge regression) or
    #' \code{"fraction"} (\eqn{\ell_1}{l1}-norm
    #' of the coefficients). Default is set to \code{"lambda"}.
    #' @param main the main title. Default is set to the model name followed
    #' by what is on the Y-axis.
    #' @param log.scale logical; indicates if a log-scale should be used
    #' when \code{xvar="lambda"}. Default is \code{TRUE}.
    #' @param standardize logical; standardize the coefficients before
    #' plotting (with the norm of the predictor). Default is \code{TRUE}.
    #' @param label vector indicating the names associated to the plotted
    #' variables. When specified, a legend is drawn in order to identify
    #' each variable. Only relevant when the number of predictor is
    #' small. Remind that the intercept does not count. Default is
    #' \code{NULL}.
    #' @param plot logical; indicates if the graph should be plotted on
    #' call. Default is \code{TRUE}.
    #'
    #' @return a \pkg{ggplot2} object which can be plotted via the
    #' \code{print} method.
    #'
    #' @examples \dontrun{
    #' ## Simulating multivariate Gaussian with blockwise correlation
    #' ## and piecewise constant vector of parameters
    #' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
    #' cor <- 0.75
    #' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
    #' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
    #' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
    #' diag(Sigma) <- 1
    #' n <- 50
    #' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
    #' y <- 10 + x %*% beta + rnorm(n,0,10)
    #'
    #' ## Plot the Lasso path
    #' plot(elastic.net(x,y, lambda2=0), main="Lasso solution path")
    #' ## Plot the Elastic-net path
    #' plot(enet, main = "Elastic-net solution path")
    #' ## Plot the Elastic-net path (fraction on X-axis, unstandardized coefficient)
    #' plot(elastic.net(x,y, lambda2=10), standardize=FALSE, xvar="fraction")
    #' ## Plot the Bounded regression path (fraction on X-axis)
    #' plot(bounded.reg(x,y, lambda2=10), xvar="fraction")
    #' }
    #'
    #' @importFrom graphics plot
    #' @import ggplot2 scales grid methods
    #' @export
    plot = function(xvar = "lambda",
                    main = paste(self$penalty," path", sep=""),
                    log.scale = TRUE, standardize=TRUE, labels = NULL, plot = TRUE, ...) {
      
      lambda <- as.numeric(unlist(self$major_penalty))
      if (length(lambda) == 1) {
        stop("Not available when the leading vector of penalties boild down to a scalar.")
      }
      
      nzeros <- which(colSums(private$beta) != 0)
      if (length(nzeros) == 0) {
        stop("Nothing to plot: all coefficients are zero.")
      }
      
      beta  <- as.matrix(private$beta[, nzeros, drop = FALSE])
      rownames(beta) <- NULL ## avoid warning message in ggplot2
      
      if (standardize) beta <- scale(beta, FALSE, 1/private$data$norm_X[nzeros])

      if (xvar == "fraction") {
        xv <-  apply(abs(beta),1,sum)/max(apply(abs(beta),1,sum))
      } else {
        xv <- lambda
      }
      
      ## Creating the data.frame fior ggploting purposes
      data.coef <- melt(data.frame(xvar=xv, beta=beta),id="xvar")
      if (is.null(labels)) {
        data.coef$labels <- factor(rep(nzeros, each=length(xv)))
      } else {
        if (sum(is.na(labels[nzeros]))>0 ) {
          labels <- NULL
          warning("The number of label is wrong, ignoring them.")
          data.coef$labels <- factor(rep(nzeros, each=length(xv)))
        } else {
          data.coef$labels <- factor(rep(labels[nzeros], each=length(xv)))
        }
      }
      colnames(data.coef) <- c("xvar","var","coef", "variables")
      d <- ggplot(data.coef,aes(x=xvar,y=coefficients, colour=variables, group=var)) +
        geom_line(aes(x=xvar,y=coef)) +  geom_hline(yintercept=0, alpha=0.5, linetype="dotted") +
        ylab(ifelse(standardize, "standardized coefficients","coefficients")) + ggtitle(main) +
        theme_bw()
      
      if (xvar=="lambda") {
        d <- d + xlab(switch(
          self$penalty,
          "ridge" = ifelse(log.scale,expression(log[10](lambda[2])),expression(lambda[2])),
                    ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1]))
          ))
        if (log.scale)
          d <- d + scale_x_log10() + annotation_logticks(sides="b")
      } else {
        d <- d + xlab(expression(paste("|",beta[lambda[1]],"|",{}[1]/max[lambda[1]],"|",beta[lambda[1]],"|",{}[1],sep="")))
      }
      
      if (is.null(labels)) {
        d <- d + theme(legend.position="none") 
      } else {
        if (length(labels[nzeros]) != length(nzeros)) {
          d <- d + theme(legend.position="none")
        }
      }
      if (plot) print(d)
      
      invisible(d)
    },
    #' Penalized criteria based on estimation of degrees of freedom
    #'
    #' Produce a plot or send back the values of some penalized criteria
    #' accompanied with the vector(s) of parameters selected
    #' accordingly. The default behavior plots the BIC and the AIC (with
    #' respective factor \eqn{\log(n)}{log(n)} and \eqn{2}{2}) yet the user can specify any
    #' penalty.
    #'
    #' @usage criteria(object, penalty=setNames(c(2, log(p)), c("AIC","BIC")), sigma=NULL,
    #'            log.scale=TRUE, xvar = "lambda", plot=TRUE)
    #'
    #' @param object output of a fitting procedure of the \pkg{quadrupen}
    #' package (e.g. \code{\link{elastic.net}}). Must be of class
    #' \code{quadrupen}.
    #' @param penalty a vector with as many penalties a desired. The
    #' default contains the penalty corresponding to the AIC and the BIC
    #' (\eqn{2}{2} and \eqn{\log(n)}{log(n)}). Setting the "names"
    #' attribute, as done in the default definition, leads to outputs
    #' which are easier to read.
    #' @param sigma scalar: an estimate of the residual variance. When
    #' available, it is plugged-in the criteria, which may be more
    #' relevant. If \code{NULL} (the default), it is estimated as usual
    #' (see details).
    #' @param xvar variable to plot on the X-axis: either \code{"df"}
    #' (the estimated degrees of freedom), \code{"lambda"}
    #' (\eqn{\lambda_1}{lambda1} penalty level) or \code{"fraction"}
    #' (\eqn{\ell_1}{l1}-norm of the coefficients). Default is set to
    #' \code{"lambda"}.
    #' @param log.scale logical; indicates if a log-scale should be used
    #' when \code{xvar="lambda"}. Default is \code{TRUE}.
    #' @param plot logical; indicates if the graph should be plotted on
    #' call. Default is \code{TRUE}.
    #'
    #' @return When \code{plot} is set to \code{TRUE}, an invisible
    #' \pkg{ggplot2} object is returned, which can be plotted via the
    #' \code{print} method. On the other hand, a list with a two data
    #' frames containing the criteria and the chosen vector of parameters
    #' are returned.
    #' @seealso \code{\linkS4class{quadrupen}}.
    #'
    #' @note When \code{sigma} is provided, the criterion takes the form
    #'
    #' \if{latex}{\deqn{\left\|\mathbf{y} - \mathbf{X} \hat{\beta} \right\|^2 +
    #' \mathrm{penalty} \times \frac{\hat{\mathrm{df}}}{n} \ \sigma^2.}}
    #' \if{html}{\out{ <center> RSS + penalty * df / n * sigma<sup>2</sup> </center>}}
    #' \if{text}{\deqn{RSS + penalty * df / n * sigma^2}}
    #'
    #' When it is unknown, it writes
    #'
    #' \if{latex}{\deqn{\log\left(\left\|\mathbf{y} - \mathbf{X} \hat{\beta} \right\|^2\right) +
    #' \mathrm{penalty} \times \hat{\mathrm{df}}.}}
    #' \if{html}{\out{ <center> n*log(RSS) + penalty * df </center>}}
    #' \if{text}{\deqn{n*log(RSS) + penalty * df}}
    #'
    #' Estimation of the degrees of freedom (for the elastic-net, the
    #' LASSO and also bounded regression) are computed by applying and
    #' adpating the results of Tibshirani and Taylor (see references
    #' below).
    #'
    #' @references Ryan Tibshirani and Jonathan Taylor. Degrees of
    #' freedom in lasso problems, Annals of Statistics, 40(2) 2012.
    #'
    #' @name criteria,quadrupen-method
    #' @aliases criteria,quadrupen-method
    #' @aliases criteria.quadrupen
    #' @aliases criteria
    #' @docType methods
    #' @rdname criteria.quadrupen
    #'
    #' @examples \dontrun{
    #' ## Simulating multivariate Gaussian with blockwise correlation
    #' ## and piecewise constant vector of parameters
    #' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
    #' cor <- 0.75
    #' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
    #' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
    #' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
    #' diag(Sigma) <- 1
    #' n <- 50
    #' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
    #' y <- 10 + x %*% beta + rnorm(n,0,10)
    #'
    #' ## Plot penalized criteria for the Elastic-net path
    #' criteria(elastic.net(x,y, lambda2=1))
    #'
    #' #' Plot penalized criteria for the Bounded regression
    #' criteria(bounded.reg(x,y, lambda2=1))
    #' }
    #'
    #' @import ggplot2 reshape2 scales grid methods
    #' @exportMethod criteria
    criteria = function(penalty=setNames(c(2, log(self$ncoef)), c("AIC","BIC")), sigma=NULL,
                         log.scale=TRUE, xvar = "lambda", plot=TRUE) {
      
      betas <- private$beta
      lambda <- as.numeric(unlist(self$major_penalty))
      
      n <- self$nsample
      p <- self$ncoef
      
      ## Compute generalized cross-validation
      GCV <- self$deviance/(n*(1 + self$degrees_freedom/n))^2
      
      ## compute the penalized criteria
      if (is.null(sigma)) {
        crit <- sapply(penalty, function(pen) log(self$deviance) + pen * self$degrees_freedom/n)
      } else {
        crit <- sapply(penalty, function(pen) self$deviance/n + pen * self$degrees_freedom/n * sigma^2)
      }
      
      ## put together all relevant information about those criteria
      criterion <- data.frame(crit, GCV=GCV, df=self$degrees_freedom, lambda=lambda, fraction = apply(abs(betas),1,sum)/max(apply(abs(betas),1,sum)), row.names=1:nrow(crit))
      
      ## recover the associated vectors of parameters
      beta.min  <- t(betas[apply(crit, 2, which.min), ])
      ##     beta.min <- cbind2(beta.min, betas[, which.min(GCV)])
      if (!is.null(dim(beta.min)))
        colnames(beta.min) <- names(penalty)
      
      ## plot the critera, if required
      if (plot) {
        
        if (length(lambda) == 1) {
          stop("Not available when the leading vector of penalties boild down to a scalar.")
        }
        
        data.plot <- melt(criterion, id=xvar, measure=1:length(penalty), variable.name="criterion", value.name="value")
        rownames(data.plot) <- 1:nrow(data.plot)
        
        colnames(data.plot)[1] <- "xvar"
        
        xlab <- switch(xvar,
            "fraction" = expression(paste("|",beta[lambda[1]],"|",{}[1]/max[lambda[1]],"|",beta[lambda[1]],"|",{}[1],sep="")),
            "df" = "Estimated degrees of freedom",
              switch(self$penalty, 
                     "ridge" = ifelse(log.scale,expression(log[10](lambda[2])),expression(lambda[2])),
                               ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1])) ) )
        
        d <- ggplot(data.plot, aes(x=xvar, y=value, colour=criterion, group=criterion)) +
          geom_line(aes(x=xvar,y=value)) + geom_point(aes(x=xvar,y=value)) +
          labs(x=xlab, y="criterion's value",  title=paste("Information Criteria for a", self$penalty,"fit"))
        
        if (log.scale & (xvar=="lambda")) {
          d <- d + scale_x_log10()
        }
        print(d)
        return(invisible(d))
      } else {
        return(list(criterion=criterion, beta.min=beta.min))
      }
      
      return(list(criterion=criterion, beta.min=beta.min))
    }
  )
)

## Auxiliary functions to check the given class of an objet
isQuadrupenFit <- function(Robject) {inherits(Robject, "QuadrupenFit"          )}

#' @export
fitted.QuadrupenFit <- function(object, ...) {
  stopifnot(isQuadrupenFit(object))
  object$fitted
}

#' @export
predict.QuadrupenFit <- function(object, newx = NULL, ...) {
  stopifnot(isQuadrupenFit(object))
  object$predict(newx = newx, ...)
}

#' @export
residuals.QuadrupenFit <- function(object, newx=NULL, newy=NULL, ...) {
  stopifnot(isQuadrupenFit(object))
  if (is.null(newx) | is.null(newy)) {
    res <- object$residuals
  } else {
    n <- length(res$lambda1)
    res <- matrix(rep(newy, n), ncol=n) - predict(object, newx)
  }
}

#' @export
deviance.QuadrupenFit <- function(object, ...) {
  stopifnot(isQuadrupenFit(object))
  object$deviance
}





#' #' @field major_penalty vector of "leading" penalties (either l1 or l2)
#' major_penalty = function(value) {
#'   switch(self$penalty,
#'          "ridge" = data.frame(lambda2=private$lambda2),
#'          data.frame(lambda1=private$lambda1))
#' },
#' #' @field major_penalty vector of "minor" penalties (either l1 or l2)
#' minor_penalty = function(value) {
#'   switch(self$penalty,
#'          "ridge" = setNames(private$lambda1, "lambda1"),
#'          setNames(private$lambda2, "lambda2"))
#' },
