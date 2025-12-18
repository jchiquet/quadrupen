
## Auxiliary functions to check the given class of an object
isQuadrupenFit <- function(Robject) {inherits(Robject, "QuadrupenFit")}

#' @export
fitted.QuadrupenFit <- function(object, ...) {
  stopifnot(isQuadrupenFit(object))
  object$fitted
}

#' @export
predict.QuadrupenFit <- function(object, newx = NULL, selection = NULL, ...) {
  stopifnot(isQuadrupenFit(object))
  object$predict(newx = newx, s = s, ...)
}

#' @export
coef.QuadrupenFit <- function(object, selection = NULL, ...) {
  stopifnot(isQuadrupenFit(object))
  object$get_model(
    selection = selection,
    type = "coefficients"
  )
}

#' @export
residuals.QuadrupenFit <- function(object, newx=NULL, newy=NULL, ...) {
  stopifnot(isQuadrupenFit(object))
  if (is.null(newx) | is.null(newy)) {
    res <- object$residuals
  } else {
    n <- length(object$major_tuning)
    res <- matrix(rep(newy, n), ncol=n) - predict(object, newx)
  }
  res
}

#' @export
deviance.QuadrupenFit <- function(object, ...) {
  stopifnot(isQuadrupenFit(object))
  object$deviance
}

criteria <- 
  function(object, 
           penalty=
             setNames(c(2, log(self$nsample), log(self$ncoef), log(self$nsample) + 2*log(self$ncoef)),
                      c("AIC","BIC", "mBIC", "eBIC")), sigma=NULL) {
    UseMethod("criteria", object)
  }

#' @export
criteria.QuadrupenFit <- 
  function(object, 
           penalty=
             setNames(c(2, log(self$nsample), log(self$ncoef), log(self$nsample) + 2*log(self$ncoef)),
                      c("AIC","BIC", "mBIC", "eBIC")), sigma=NULL) {
    stopifnot(isQuadrupenFit(object))
    object$criteria(penalty, sigma)
    object$information_criteria
}

#' @export
cross_validate <- 
  function(object, 
           K       = 10,
           folds   = split(sample(1:self$nsample), rep(1:K, length=self$nsample)),
           lambda2 = self$minor_tuning, verbose = TRUE, cores = parallel::detectCores() - 2
  ) {
    UseMethod("cross_validate", object)
  }

#' @export
cross_validate.QuadrupenFit <- 
  function(object, 
           K       = 10,
           folds   = split(sample(1:self$nsample), rep(1:K, length=self$nsample)),
           lambda2 = self$minor_tuning, verbose = TRUE, cores = parallel::detectCores() - 2
    ) {
      stopifnot(isQuadrupenFit(object))
  object$cross_validate(K, folds, lambda2, verbose, cores)
  object$cross_validation
}

#' @export
stability <- 
  function(object, n_subsamples   = 50,
           subsample_size = floor(self$nsample/2),
           subsamples     = replicate(n_subsamples, sample(1:self$nsample, subsample_size), simplify=FALSE),
           weakness       = 1,
           verbose        = TRUE,
           cores          = parallel::detectCores() - 2) {
    UseMethod("stability", object)
  }

#' @export
stability.QuadrupenFit <- 
  function(object, n_subsamples   = 50,
           subsample_size = floor(self$nsample/2),
           subsamples     = replicate(n_subsamples, sample(1:self$nsample, subsample_size), simplify=FALSE),
           weakness       = 1,
           verbose        = TRUE,
           cores          = parallel::detectCores() - 2) {
    stopifnot(isQuadrupenFit(object))
  object$stability(n_subsamples, subsample_size, subsamples, weakness, verbose, cores)
  object$stability_path
}
