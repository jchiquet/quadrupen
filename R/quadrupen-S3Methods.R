
#' Auxiliary functions to check the given class of an object
#' @param Robject an R object to evaluate
#' @return logical
isQuadrupenFit <- function(Robject) {inherits(Robject, "QuadrupenFit")}

#' Extracts model fitted values
#'
#' @param object a [`QuadrupenFit`] object
#' @param ... not used, only here for S3 compatibility
#'
#' @return A matrix of fitted values extracted from `object`.
#'
#' @export
fitted.QuadrupenFit <- function(object, ...) {
  stopifnot(isQuadrupenFit(object))
  object$fitted
}

#' Perform model prediction 
#' 
#' @description Predict response for new sample based on the current model
#' 
#' @param object an R object to evaluate
#' @param newx matrix of new values for the regressor with which to predict. If omitted, the fitted values are used.
#' @param selection either a character (model selection criteria) of a scalar (lambda value)
#' @param ... not used, only here for S3 compatibility
#' 
#' @return a vector of predicted value
#' @export
predict.QuadrupenFit <- function(object, newx = NULL, selection = NULL, ...) {
  stopifnot(isQuadrupenFit(object))
  object$predict(newx = newx, selection = selection)
}

#' Extract model coefficients
#'
#' @description Extracts model coefficients from a [`QuadrupenFit`] object
#' @param object a [`QuadrupenFit`] object
#' @param selection either a character (model selection criteria) of a scalar (lambda value)
#' @param ... not used, only here for S3 compatibility
#' @return a vector of coefficients
#' @export
coef.QuadrupenFit <- function(object, selection = NULL, ...) {
  stopifnot(isQuadrupenFit(object))
  object$get_model(
    selection = selection,
    type = "coefficients"
  )
}

#' Extract model residuals 
#' 
#' @description Extracts model residuals from a [`QuadrupenFit`] object
#' 
#' @param object a [`QuadrupenFit`] object
#' @param newx matrix of new values for the regressor with which to predict. If omitted, the fitted values are used.
#' @param newy vector of new values for the response with which to compute the residuals. If omitted, the fitted values are used.
#' @param ... not used, only here for S3 compatibility
#' @return Matrix of residuals, each column corresponding to a value of \code{lambda1}.
#' @export
residuals.QuadrupenFit <- function(object, newx=NULL, newy, ...) {
  stopifnot(isQuadrupenFit(object))
  if (is.null(newx) | is.null(newy)) {
    res <- object$residuals
  } else {
    n <- length(object$major_tuning)
    res <- matrix(rep(newy, n), ncol=n) - predict(object, newx)
  }
  res
}

#' Extract model deviance
#' 
#' @description Extracts the deviance of a [`QuadrupenFit`] object
#' 
#' @param object a [`QuadrupenFit`] object
#' @param ... not used, only here for S3 compatibility
#' @return a scalar
#' @export
deviance.QuadrupenFit <- function(object, ...){
  stopifnot(isQuadrupenFit(object))
  object$deviance
}

#' Penalized criteria based on estimation of degrees of freedom
#'
#' @description Produce a plot or send back the values of some penalized criteria
#' accompanied with the vector(s) of parameters selected
#' accordingly. The default behavior plots the BIC and the AIC (with
#' respective factor \eqn{\log(n)}{log(n)} and \eqn{2}{2}) yet the user can specify any
#' penalty.
#'
#' @param object output of a fitting procedure of the \pkg{quadrupen}
#' package (e.g. [elastic.net()]).
#' @param penalty a vector with as many penalties a desired. The
#' default contains the penalty corresponding to the AIC and the BIC
#' (\eqn{2}{2} and \eqn{\log(n)}{log(n)}). Setting the "names"
#' attribute, as done in the default definition, leads to outputs
#' which are easier to read.
#' @param sigma scalar: an estimate of the residual variance. When
#' available, it is plugged-in the criteria, which may be more
#' relevant. If `NULL` (the default), it is estimated as usual
#' (see details).
#' 
#' @return an object with class [`Criteria`] is sent back and stored as a 
#' field of the original [`QuadrupenFit`] object.
#'
#' @note When `sigma` is provided, the criterion takes the form
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
#' @importFrom stats setNames
#' @export
criteria <- 
  function(object, 
           penalty=
             setNames(c(2, log(object$nsample), log(object$ncoef), log(object$nsample) + 2*log(object$ncoef)),
                      c("AIC","BIC", "mBIC", "eBIC")), sigma=NULL) {
    UseMethod("criteria", object)
  }

#' @describeIn criteria S3 method for information criteria  of a [`QuadrupenFit`]
#' @export
criteria.QuadrupenFit <- 
  function(object, 
           penalty=
             setNames(c(2, log(object$nsample), log(object$ncoef), log(object$nsample) + 2*log(object$ncoef)),
                      c("AIC","BIC", "mBIC", "eBIC")), sigma=NULL) {
    stopifnot(isQuadrupenFit(object))
    object$criteria(penalty, sigma)
    object$information_criteria
}

#' Cross-validation for Quadrupen object
#' 
#' @description Function that computes K-fold cross-validated error of a
#' \code{quadrupen} fit, possibly on a grid of `lambda1`, `lambda2`.
#'
#' @param object an R6 object with class [`QuadrupenFit`]
#'  
#' @param K integer indicating the number of folds. Default is 10.
#'
#' @param folds list of `K` vectors that describes the folds to
#' use for the cross-validation. By default, the folds are randomly
#' sampled with the specified K. The same folds are used for each
#' values of `lambda2`.
#'
#' @param lambda2 tunes the \eqn{\ell_2}{l2}-penalty (ridge-like) of
#' the fit. If none is provided, a vector of values is generated and
#' a CV is performed on a grid of `lambda2` and `lambda1`,
#' using the same folds for each `lambda2`.
#'
#' @param verbose logical; indicates if the progression (the current
#' `lambda2` should be displayed. Default is `TRUE.`
#'
#' @param cores the number of cores to use. The default uses all
#' the cores available.
#'
#' @note If the user runs the fitting method with option
#' \code{'bulletproof'} set to \code{FALSE}, the algorithm may stop
#' at an early stage of the path. Early stops are handled internally,
#' in order to provide results on the same grid of penalty tuned by
#' \eqn{\lambda_1}{lambda1}.  This is done by means of \code{NA}
#' values, so as mean and standard error are consistently
#' evaluated. If, while cross-validating, the procedure experiences
#' too much early stoppings, a warning is sent to the user, in which
#' case you should reconsider the grid of \code{lambda1} used for the
#' cross-validation.  If \code{bulletproof} is \code{TRUE} (the
#' default), there is nothing to worry about, except a possible slow
#' down when any switching to the proximal algorithm is required.
#'
#' @return an object with class [`CrossValidation`] is sent back and stored as a 
#' field of the original [`QuadrupenFit`] object.
#'
#' @examples \dontrun{
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
#' cor  <- 0.75
#' Soo  <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variable
#' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
#' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.1
#' diag(Sigma) <- 1
#' n <- 100
#' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
#' y <- 10 + x %*% beta + rnorm(n,0,10)
#' 
#' enet <- elastic.net(x, y, nlambda1=50)
#' 
#' ## Use fewer lambda1 values by overwritting the default parameters
#' ## and cross-validate over the sequences lambda1 and lambda2
#' cv.grid <- cross_validate(enet, lambda2=10^seq(2,-2,len=50))
#' ## Rerun simple cross-validation with the appropriate lambda2
#' cv.10K <- crossval(x,y, lambda2=cv.grid$lambda2_min)
#' ## Try leave one out also
#' cv.loo <- crossval(x,y, K=n, lambda2=cv.grid$lambda2_min)
#'
#' plot(cv.grid)
#' plot(cv.10K)
#' plot(cv.loo)
#'
#' ## Performance for selection purpose
#' cat("\nFalse positives with the minimal 10-CV choice: ", sum(sign(beta) != sign(cv.10K$beta_min )))
#' cat("\nFalse positives with the minimal LOO-CV choice: ", sum(sign(beta) != sign(cv.loo$beta_min)))
#' }
#' 
#' @export
cross_validate <- 
  function(object, 
           K       = 10,
           folds   = split(sample(1:object$nsample), rep(1:K, length=object$nsample)),
           lambda2 = object$minor_tuning, verbose = TRUE, cores = parallel::detectCores() - 2
  ) {
    UseMethod("cross_validate", object)
  }

#' @describeIn cross_validate S3 method for cross-validation  of a [`QuadrupenFit`]
#' @export
cross_validate.QuadrupenFit <- 
  function(object, 
           K       = 10,
           folds   = split(sample(1:object$nsample), rep(1:K, length=object$nsample)),
           lambda2 = object$minor_tuning, verbose = TRUE, cores = parallel::detectCores() - 2
    ) {
      stopifnot(isQuadrupenFit(object))
  object$cross_validate(K, folds, lambda2, verbose, cores)
  object$cross_validation
}

#' Stability selection for Quadrupen object
#' 
#' @description Compute the stability path of a (possibly randomized) fitting
#' procedure as introduced by Meinshausen and Buhlmann (2010).
#'
#' @param object an R6 object with class [`QuadrupenFit`]
#' 
#' @param n_subsamples integer indicating the number of subsamplings
#' used to estimate the selection probabilities. Default is 100.
#'
#' @param subsample_size integer indicating the size of each subsamples.
#' Default is `floor(n/2)`.
#'
#' @param subsamples list with `subsamples` entries with vectors
#' describing the folds to use for the stability procedure. By
#' default, the folds are randomly sampled with the specified
#' \code{n_subsamples} and `subsample_size` argument.
#'
#' @param weakness Coefficient used for randomizing the weights of each features.
#' Default is 1` for no randomization. See details below.
#' 
#' @param verbose logical; indicates if the progression should be
#' displayed. Default is `TRUE`.
#'
#' @param cores the number of cores to use. The default uses all
#' the cores available.
#'
#' @return an object with class [`StabilityPath`] is sent back and stored as a 
#' field of the original [`QuadrupenFit`] object.
#'
#' @note When `weakness < 1`, the \code{penscale} argument
#' that weights the penalty tuned by \eqn{\lambda_1}{lambda1} is
#' perturbed (divided) for each subsample by a random variable
#' uniformly distributed on
#' \if{latex}{\eqn{[\alpha,1]}}\if{html}{[&#945;,1]}\if{text}{\eqn{[alpha,1]}},
#' where
#' \if{latex}{\eqn{\alpha}}\if{html}{&#945;}\if{text}{\eqn{alpha}} is
#' the weakness parameter.
#'
#' If the user runs the fitting method with option
#' `'bulletproof'` set to `FALSE`, the algorithm may stop
#' at an early stage of the path. Early stops of the underlying
#' fitting function are handled internally, in the following way: we
#' chose to simply skip the results associated with such runs, in
#' order not to bias the stability selection procedure. If it occurs
#' too often, a warning is sent to the user, in which case you should
#' reconsider the grid of `lambda1` for stability selection. If
#' `bulletproof` is `TRUE` (the default), there is nothing
#' to worry about, except a possible slow down when any switching to
#' the proximal algorithm is required.
#'
#' @references N. Meinshausen and P. Buhlmann (2010). Stability
#' Selection, JRSS(B).
#'
#' @examples \dontrun{
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
#' Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
#' Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
#' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
#' diag(Sigma) <- 1
#' n <- 100
#' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
#' y <- 10 + x %*% beta + rnorm(n,0,10)
#'
#' ## Build a vector of label for true nonzeros
#' labels <- rep("irrelevant", length(beta))
#' labels[beta != 0] <- c("relevant")
#' labels <- factor(labels, ordered=TRUE, levels=c("relevant","irrelevant"))
#'
#' ## Call to stability selection function, 200 subsampling
#' stab <- stability(x,y, subsamples=200, lambda2=1, min.ratio=1e-2)
#' ## Recover the selected variables for a given cutoff
#' ## and per-family error rate, without producing any plot
#' stabpath <- plot(stab, cutoff=0.75, PFER=1, plot=FALSE)
#'
#' cat("\nFalse positives for the randomized Elastic-net with stability selection: ",
#'      sum(labels[stabpath$selected] != "relevant"))
#' cat("\nDONE.\n")
#'}
#'
#' @export
stability <- 
  function(object, n_subsamples   = 50,
           subsample_size = floor(object$nsample/2),
           subsamples     = replicate(n_subsamples, sample(1:object$nsample, subsample_size), simplify=FALSE),
           weakness       = 1,
           verbose        = TRUE,
           cores          = parallel::detectCores() - 2) {
    UseMethod("stability", object)
  }

#' @describeIn stability S3 method for stability selection of a [`QuadrupenFit`]
#' @export
stability.QuadrupenFit <- 
  function(object, n_subsamples   = 50,
           subsample_size = floor(object$nsample/2),
           subsamples     = replicate(n_subsamples, sample(1:object$nsample, subsample_size), simplify=FALSE),
           weakness       = 1,
           verbose        = TRUE,
           cores          = parallel::detectCores() - 2) {
    stopifnot(isQuadrupenFit(object))
  object$stability(n_subsamples, subsample_size, subsamples, weakness, verbose, cores)
  object$stability_path
}

#' Plot method
#' 
#' @description Produce a plot of the solution path of a \code{quadrupen} fit.
#' 
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
