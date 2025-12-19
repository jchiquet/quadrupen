#' Class Criteria
#' 
#' Class of object returned by a call to the \code{$criteria()} method.
#' 
#' \@param lambda vector of \eqn{\lambda_1}{lambda}
#' (\eqn{\ell_1}{l1} or \eqn{\ell_\infty}{infinity} penalty levels)
#' for which each cross-validation has been performed.
#' @param lambda2 vector (or scalar) of \eqn{\ell_2}{l2}-penalty levels for
#' which each cross-validation has been performed.
#' 
#' @importFrom dplyr filter
#' 
#' @export
Criteria <- R6::R6Class(
  classname = "CrossValidation",
  private = list(
    value     = NA
  ),
  active = list(
    lambda = function(value) {private$value$lambda},
    data   = function(value) {private$value},
    names  = function(value) {colnames(private$value)[1:(ncol(private$value)-3)]}
  ),
  public = list(
    initialize = function(value) {
      private$value <- value
    },
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
    plot = function(log.scale=TRUE, xvar = c("lambda", "fraction", "df"), title = "Information Criteria") {
      if (length(self$lambda) == 1) {
        stop("Not available when the leading vector of penalties boild down to a scalar.")
      }
      xvar <- match.arg(xvar)

      data.plot <- melt(self$data, id=xvar, measure=1:(ncol(self$data)-3), variable.name="criterion", value.name="value")
      rownames(data.plot) <- 1:nrow(data.plot)
      colnames(data.plot)[1] <- "xvar"
      
      
      data_plot <- 
        dplyr::select(self$data, (1:(ncol(self$data) - 3)) | starts_with(xvar) ) |> 
                    tidyr::pivot_longer(cols = -xvar, names_to = "criterion") |>
                    dplyr::rename(xvar = 1)

      xlab <- switch(
        xvar,
        "fraction" = expression(paste("|",beta[lambda[1]],"|",{}[1]/max[lambda[1]],"|",beta[lambda[1]],"|",{}[1],sep="")),
        "df" = "Estimated degrees of freedom",
        ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1]))
        )
      
      d <- ggplot(data_plot) +
        aes(x=xvar, y=value, colour=criterion, group=criterion) +
        geom_line() + geom_point() + theme_bw() + 
        labs(x = xlab, y = "criterion's value",  title = title)
      
      if (log.scale & (xvar == "lambda")) {
        d <- d + scale_x_log10()
      }
      print(d)
      return(invisible(d))
    }
  )
)  
  