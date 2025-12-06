
## Auxiliary functions to check the given class of an object
isQuadrupenFit <- function(Robject) {inherits(Robject, "QuadrupenFit")}

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
