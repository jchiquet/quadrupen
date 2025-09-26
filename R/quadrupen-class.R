#' Class "quadrupen"
#'
#' Class of object returned by any fitting function of the
#' \pkg{quadrupen} package (\code{elastic.net} or
#' \code{bounded.reg}).
#'
#' @section Slots: \describe{
#'
#' \item{\code{coefficients}:}{Matrix (class \code{"dgCMatrix"}) of
#' coefficients with respect to the original input. The number of
#' rows corresponds the length of \code{lambda1}.}
#'
#' \item{\code{active.set}:}{Matrix (class \code{"dgCMatrix"}, generally
#' sparse) indicating the 'active' variables, in the sense that they
#' activate the constraints. For the \code{\link{elastic.net}}, it
#' corresponds to the nonzero variables; for the
#' \code{\link{bounded.reg}} function, it is the set of variables
#' reaching the boudary along the path of solutions.}
#'
#' \item{\code{intercept}:}{logical; indicates if an intercept has
#'  been included to the model.}
#'
#' \item{\code{mu}:}{A vector (class \code{"numeric"})
#' containing the successive values of the (unpenalized) intercept.
#' Equals to zero if \code{intercept} has been set to \code{FALSE}.}
#'
#' \item{\code{normx}:}{Vector (class \code{"numeric"}) containing the
#' square root of the sum of squares of each column of the design
#' matrix.}
#'
#' \item{\code{penscale}:}{Vector \code{"numeric"} with real positive
#' values that have been used to weight the penalty tuned by
#' \eqn{\lambda_1}{lambda1}.}
#'
#' \item{\code{penalty}:}{Object of class \code{"character"}
#' indicating the method used (\code{"elastic-net"} or \code{"bounded
#' regression"}).}
#'
#' \item{\code{naive}:}{logical; was the \code{naive} mode on?}
#'
#' \item{\code{lambda1}:}{Vector (class \code{"numeric"}) of penalty
#' levels (either \eqn{\ell_1}{l1} or \eqn{\ell_\infty}{l-infinity})
#' for which the model has eventually been fitted.}
#'
#' \item{\code{lambda2}:}{Scalar (class \code{"numeric"}) for the
#' amount of \eqn{\ell_2}{l2} (ridge-like) penalty.}
#'
#' \item{\code{control}:}{Object of class \code{"list"} with low
#' level options used for optimization.}
#'
#' \item{\code{monitoring}:}{List (class \code{"list"}) which
#' contains various indicators dealing with the optimization
#' process.}
#'
#' \item{\code{residuals}:}{Matrix of residuals, each column
#' corresponding to a value of \code{lambda1}.}
#'
#' \item{\code{df}:}{Estimated degree of freedoms for the successive
#' \code{lambda1}.  Only available for 'elastic.net' using tCholesky
#' factorization.}
#'
#' \item{\code{r.squared}:}{Vector (class \code{"numeric"}) given the
#' coefficient of determination as a function of lambda1.}
#'
#' \item{\code{fitted}:}{Matrix of fitted values, each column
#' corresponding to a value of \code{lambda1}.}  }
#'
#' @section Methods:
#' This class comes with the usual \code{predict(object, newx, ...)},
#' \code{fitted(object, ...)}, \code{residuals(object, ...)},
#' \code{print(object, ...)}, \code{show(object)} and
#' \code{deviance(object, ...)} generic (undocumented) methods.
#'
#' A specific plotting method is available and documented
#' (\code{\link{plot.quadrupen}}).
#'
#' @aliases fitted,quadrupen-method predict,quadrupen-method
#' deviance,quadrupen-method print,quadrupen-method
#' show,quadrupen-method residuals,quadrupen-method
#'
#' @docType class
#'
#' @keywords class
#'
#' @seealso See also \code{\link{plot.quadrupen}}.
#'
#' @name quadrupen-class
#' @rdname quadrupen-class
#'
#' @exportClass quadrupen
#' @exportMethod fitted
#' @exportMethod residuals
#' @exportMethod predict
#' @exportMethod deviance
#' @exportMethod print
#' @exportMethod show
#'
#' @importFrom stats fitted predict residuals deviance
#'
setClassUnion("mat", c("Matrix","matrix"))
setClassUnion("naive", c("NULL","logical"))
setClass("quadrupen",
  representation = representation(
     coefficients  = "Matrix",
     active.set    = "Matrix",
     intercept     = "logical"  ,
     mu            = "numeric"  ,
     normx         = "numeric"  ,
     fitted        = "mat"      ,
     residuals     = "mat"      ,
     df            = "numeric"  ,
     r.squared     = "numeric"  ,
     penscale      = "numeric"  ,
     penalty       = "character",
     naive         = "naive"    ,
     lambda1       = "numeric"  ,
     lambda2       = "numeric"  ,
     control       = "list"     ,
     monitoring    = "list")
)

setMethod("print", "quadrupen", definition =
   function(x, ...) {
     ncoef <- ncol(x@coefficients)
     if (!is.null(x@naive)) {
       if (x@naive) {
         cat("Linear regression with", x@penalty, "penalizer, no rescaling of the coefficients (naive).\n")
       } else {
         cat("Linear regression with", x@penalty, "penalizer, coefficients rescaled by (1+lambda2).\n")
       }
     } else {
       cat("Linear regression with", x@penalty, "penalizer.\n")
     }
     if (x@intercept) {
       cat("- number of coefficients:", ncoef,"+ intercept\n")
     } else {
       cat("- number of coefficients:", ncoef,"(no intercept)\n")
     }

     major <- major.pen(x)
     minor <- minor.pen(x)
     cat("- penalty parameter ",names(major), ": ", length(major), " points from ",
         format(max(major), digits = 3)," to ",
         format(min(major), digits = 3),"\n", sep="")
     cat("- penalty parameter ",names(minor),": ", minor, "\n", sep="")

     invisible(x)
   }
)

## this is simply an internal basic function to get vector if "leading" penalties (either l1 or l2)
setGeneric("major.pen", function(object) {standardGeneric("major.pen")})
setMethod("major.pen", "quadrupen", definition =
   function(object) {
     switch(object@penalty, "ridge" = data.frame(lambda2=object@lambda2), data.frame(lambda1=object@lambda1))
   })
## this is simply an internal basic function to get vector if "minor" penalties (either l1 or l2)
setGeneric("minor.pen", function(object) {standardGeneric("minor.pen")})
setMethod("minor.pen", "quadrupen", definition =
   function(object) {
     switch(object@penalty, "ridge" = setNames(object@lambda1, "lambda1"), setNames(object@lambda2, "lambda2"))
   })

setMethod("show", "quadrupen", definition =
   function(object) {print(object)}
)

setMethod("fitted", "quadrupen", definition =
   function(object, ...) {
     return(object@fitted)
   }
)

setMethod("predict", "quadrupen", definition =
   function (object, newx=NULL, ...)  {
     if (is.null(newx)) {
       return(object@fitted)
     } else {
       return(sweep(newx %*% t(object@coefficients),2L,-object@mu,check.margin=FALSE))
     }
   }
)

setMethod("residuals", "quadrupen", definition =
   function(object, newx=NULL, newy=NULL, ...) {
     if (is.null(newx) | is.null(newy)) {
       return(object@residuals)
     } else {
       n <- length(object@lambda1)
       return(matrix(rep(newy, n), ncol=n) - predict(object, newx))
     }
   }
)

setMethod("deviance", "quadrupen", definition =
   function(object, newx=NULL, newy=NULL, ...) {
     dev <- colSums(residuals(object, newx, newy)^2)
     return(dev)
   }
)

#' Plot method for a quadrupen object
#'
#' Produce a plot of the solution path of a \code{quadrupen} fit.
#'
#' @usage plot.quadrupen(x, y, xvar = "lambda",
#'         main = paste(slot(x, "penalty")," path", sep=""),
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
#' @seealso \code{\linkS4class{quadrupen}}.
#'
#' @name plot,quadrupen-method
#' @aliases plot,quadrupen-method
#' @aliases plot.quadrupen
#' @docType methods
#' @rdname plot.quadrupen
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
#' @exportMethod plot
#' @import ggplot2 scales grid methods
#' @export
setMethod("plot", "quadrupen", definition =
   function(x, y, xvar = "lambda",
            main = paste(slot(x, "penalty")," path", sep=""),
            log.scale = TRUE, standardize=TRUE, labels = NULL, plot = TRUE, ...) {

     lambda <- as.numeric(unlist(major.pen(x)))
     if (length(lambda) == 1) {
       stop("Not available when the leading vector of penalties boild down to a scalar.")
     }

     nzeros <- which(colSums(x@coefficients) != 0)
     if (length(nzeros) == 0) {
       stop("Nothing to plot: all coefficients are zero.")
     }

     beta  <- as.matrix(x@coefficients[, nzeros])
     rownames(beta) <- NULL ## avoid warning message in ggplot2

     if (standardize) {
       beta  <- scale(beta, FALSE, 1/x@normx[nzeros])
     }

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
         ylab(ifelse(standardize, "standardized coefficients","coefficients")) + ggtitle(main)

     if (xvar=="lambda") {
       d <- d + xlab(switch(x@penalty, "ridge" = ifelse(log.scale,expression(log[10](lambda[2])),expression(lambda[2])),
                            ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1]))))
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

   }
)

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
setGeneric("criteria", function(object, penalty=setNames(c(2, log(ncol(object@coefficients))), c("AIC","BIC")), sigma=NULL,
            log.scale=TRUE, xvar = "lambda", plot=TRUE)
           {standardGeneric("criteria")})

setMethod("criteria", "quadrupen", definition =
   function(object, penalty=setNames(c(2, log(ncol(object@coefficients))), c("AIC","BIC")), sigma=NULL,
            log.scale=TRUE, xvar = "lambda", plot=TRUE) {

     betas <- object@coefficients
     lambda <- as.numeric(unlist(major.pen(object)))

     n <- nrow(residuals(object))
     p <- ncol(betas)

     ## Compute generalized cross-validation
     GCV <- deviance(object)/(n*(1+object@df/n))^2

     ## compute the penalized criteria
     if (is.null(sigma)) {
       crit <- sapply(penalty, function(pen) log(deviance(object)) + pen * object@df/n)
     } else {
       crit <- sapply(penalty, function(pen) deviance(object)/n + pen * object@df/ n * sigma^2)
     }


     ## put together all relevant information about those criteria
     criterion <- data.frame(crit, GCV=GCV, df=object@df, lambda=lambda, fraction = apply(abs(betas),1,sum)/max(apply(abs(betas),1,sum)), row.names=1:nrow(crit))

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
                      switch(object@penalty, "ridge" = ifelse(log.scale,expression(log[10](lambda[2])),expression(lambda[2])),
                             ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1])) ) )

       d <- ggplot(data.plot, aes(x=xvar, y=value, colour=criterion, group=criterion)) +
         geom_line(aes(x=xvar,y=value)) + geom_point(aes(x=xvar,y=value)) +
           labs(x=xlab, y="criterion's value",  title=paste("Information Criteria for a", slot(object, "penalty"),"fit"))

       if (log.scale & xvar=="lambda") {
         d <- d + scale_x_log10()
       }
       print(d)
       return(invisible(d))
     } else {
       return(list(criterion=criterion, beta.min=beta.min))
     }

     return(list(criterion=criterion, beta.min=beta.min))
   })

