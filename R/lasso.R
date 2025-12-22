#' Fit a linear model with lasso regularization
#'
#' Adjust a linear model with lasso regularization, that is a
#' (possibly weighted) \eqn{\ell_1}{l1}-norm. The solution path is
#' computed at a grid of values for the \eqn{\ell_1}{l1}-penalty. See
#' details for the criterion optimized.
#'
#' @inheritParams elastic.net
#'
#' @return an object with class [`QuadrupenFit`].
#'
#' @note The optimized criterion is the following: \if{latex}{\deqn{%
#' \hat{\beta}_{\lambda_1} = \arg \min_{\beta} \frac{1}{2} (y - X
#' \beta)^T (y - X \beta) + \lambda_1 \|D \beta \|_{1}, }}
#' \if{html}{\out{ <center> &beta;<sup>hat</sup>
#' <sub>&lambda;<sub>1</sub></sub> =
#' argmin<sub>&beta;</sub> 1/2 RSS(&beta) + &lambda;<sub>1</sub>
#' &#124; D &beta; &#124;<sub>1</sub>, </center> }}
#' \if{text}{\deqn{beta.hat(lambda1) = argmin_beta 1/2
#' RSS(beta) + lambda1 |D beta|1,}} where
#' \eqn{D}{D} is a diagonal matrix, whose diagonal terms are provided
#' as a vector by the \code{penscale} argument.
#'
#' @seealso See also [`QuadrupenFit`],
#' 
#' @keywords models, regression
#'
#' @examples
#' ## Simulating multivariate Gaussian with blockwise correlation
#' ## and piecewise constant vector of parameters
#' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
#' cor <- 0.75
#' Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
#' Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
#' Sigma <- Matrix::bdiag(Soo,Sww,Soo,Sww,Soo)
#' diag(Sigma) <- 1
#' n <- 50
#' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
#' y <- 10 + x %*% beta + rnorm(n,0,10)
#'
#' labels <- rep("irrelevant", length(beta))
#' labels[beta != 0] <- "relevant"
#' ## The solution path of the LASSO
#' plot(lasso(x,y), label=labels)
#'
#' @export
lasso <- function(x,
                  y,
                  lambda1   = NULL,
                  penscale  = rep(1,ncol(x)),
                  intercept = TRUE,
                  normalize = TRUE,
                  nlambda1  = ifelse(is.null(lambda1),100,length(lambda1)),
                  min.ratio = ifelse(nrow(x) <= ncol(x), 1e-2, 1e-4),
                  max.feat  = min(nrow(x),ncol(x)),
                  beta0     = numeric(ncol(x)),
                  control   = list()) {
  
  out <- elastic.net(x,
                     y,
                     lambda1   = lambda1,
                     lambda2   = 0,
                     penscale  = penscale,
                     intercept = intercept,
                     normalize = normalize,
                     nlambda1  = nlambda1,
                     min.ratio = min.ratio,
                     max.feat  = max.feat,
                     beta0     = beta0,
                     control   = control)
  # out@penalty <- "lasso"
  out

}

