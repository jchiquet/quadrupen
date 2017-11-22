\docType{methods}
\name{criteria,quadrupen-method}
\alias{criteria}
\alias{criteria.quadrupen}
\alias{criteria,quadrupen-method}
\title{Penalized criteria based on estimation of degrees of freedom}
\arguments{
  \item{object}{output of a fitting procedure of the
  \pkg{quadrupen} package (e.g. \code{\link{elastic.net}}).
  Must be of class \code{quadrupen}.}

  \item{penalty}{a vector with as many penalties a desired.
  The default contains the penalty corresponding to the AIC
  and the BIC (\eqn{2}{2} and \eqn{\log(n)}{log(n)}).
  Setting the "names" attribute, as done in the default
  definition, leads to outputs which are easier to read.}

  \item{sigma}{scalar: an estimate of the residual
  variance. When available, it is plugged-in the criteria,
  which may be more relevant. If \code{NULL} (the default),
  it is estimated as usual (see details).}

  \item{xvar}{variable to plot on the X-axis: either
  \code{"df"} (the estimated degrees of freedom),
  \code{"lambda"} (\eqn{\lambda_1}{lambda1} penalty level)
  or \code{"fraction"} (\eqn{\ell_1}{l1}-norm of the
  coefficients). Default is set to \code{"lambda"}.}

  \item{log.scale}{logical; indicates if a log-scale should
  be used when \code{xvar="lambda"}. Default is
  \code{TRUE}.}

  \item{plot}{logical; indicates if the graph should be
  plotted on call. Default is \code{TRUE}.}
}
\value{
  When \code{plot} is set to \code{TRUE}, an invisible
  \pkg{ggplot2} object is returned, which can be plotted
  via the \code{print} method. On the other hand, a list
  with a two data frames containing the criteria and the
  chosen vector of parameters are returned.
}
\description{
  Produce a plot or send back the values of some penalized
  criteria accompanied with the vector(s) of parameters
  selected accordingly. The default behavior plots the BIC
  and the AIC (with respective factor \eqn{\log(n)}{log(n)}
  and \eqn{2}{2}) yet the user can specify any penalty.
}
\note{
  When \code{sigma} is provided, the criterion takes the
  form

  \if{latex}{\deqn{\left\|\mathbf{y} - \mathbf{X}
  \hat{\beta} \right\|^2 + \mathrm{penalty} \times
  \frac{\hat{\mathrm{df}}}{n} \ \sigma^2.}} \if{html}{\out{
  <center> RSS + penalty * df / n * sigma<sup>2</sup>
  </center>}} \if{text}{\deqn{RSS + penalty * df / n *
  sigma^2}}

  When it is unknown, it writes

  \if{latex}{\deqn{\log\left(\left\|\mathbf{y} - \mathbf{X}
  \hat{\beta} \right\|^2\right) + \mathrm{penalty} \times
  \hat{\mathrm{df}}.}} \if{html}{\out{ <center> n*log(RSS)
  + penalty * df </center>}} \if{text}{\deqn{n*log(RSS) +
  penalty * df}}

  Estimation of the degrees of freedom (for the
  elastic-net, the LASSO and also bounded regression) are
  computed by applying and adpating the results of
  Tibshirani and Taylor (see references below).
}
\examples{
\dontrun{
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
cor <- 0.75
Soo <- toeplitz(cor^(0:(25-1))) ## Toeplitz correlation for irrelevant variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) \%*\% chol(Sigma))
y <- 10 + x \%*\% beta + rnorm(n,0,10)

## Plot penalized criteria for the Elastic-net path
criteria(elastic.net(x,y, lambda2=1))

##' Plot penalized criteria for the Bounded regression
criteria(bounded.reg(x,y, lambda2=1))
}
}
\references{
  Ryan Tibshirani and Jonathan Taylor. Degrees of freedom
  in lasso problems, Annals of Statistics, 40(2) 2012.
}
\seealso{
  \code{\linkS4class{quadrupen}}.
}

