#' Class Criteria
#' 
#' Class of object returned by the [`QuadrupenFit$criteria()`] method or the
#' [`criteria()`] function.  Owns [print()] and [plot()] methods.
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
    #' @field data a data frame containing the values of various information
    #' criteria (AIC, BIC, lmBIC, eBIC, GCV) along the value of \code{lambda1}.
    data   = function(value) {private$value},
    #' @field lambda vector of \eqn{\lambda_1}{lambda}
    #' (\eqn{\ell_1}{l1} or \eqn{\ell_\infty}{infinity} penalty levels)
    #' for which each cross-validation has been performed.
    lambda = function(value) {private$value$lambda},
    #' @field names a vector of characters storing the names of the precomputed criteria
    names  = function(value) {colnames(private$value)[1:(ncol(private$value) - 3)]}
  ),
  public = list(
    #' @description Constructor for a [`Criteria`] object
    #' Should be called internally by an object [`QuadrupenFit$criteria()`]
    #' @param value data frame storing output of [`QuadrupenFit$criteria()`]
    initialize = function(value) {
      private$value <- value
    },
    #' @description User friendly print method
    show = function() {
      cat("Information criteria for a Quadrupen fit.\n")
      cat("- Criteria considered: ", self$names, "\n")
      cat("- main penalty parameter lambda:", length(self$lambda), "points from",
          format(max(self$lambda), digits = 3),"to",
          format(min(self$lambda), digits = 3),"\n")
      invisible(self)
    },
    #' @description Plot the the desired criteria 
    #' 
    #' @param criteria a vector of character with the criteria to plot. 
    #' The default plot all the criteria available (stored in the field `names`)
    #' 
    #' @param log_scale logical; indicates if a log-scale should be used
    #' when `xvar="lambda"`. Default is `TRUE`.
    #' 
    #' @param xvar variable to plot on the X-axis: either `"df"`
    #' (the estimated degrees of freedom), `"lambda"`
    #' (\eqn{\lambda_1}{lambda1} penalty level) or `"fraction"`
    #' (\eqn{\ell_1}{l1}-norm of the coefficients). Default is set to
    #' `"lambda"`.
    #' 
    #' @param title graph title
    #'
    #' @return a \pkg{ggplot2} object
    plot = function(log_scale=TRUE, criteria = self$names, xvar = c("lambda", "fraction", "df"), title = "Information Criteria") {
      if (length(self$lambda) == 1) {
        stop("Not available when the leading vector of penalties boild down to a scalar.")
      }
      xvar <- match.arg(xvar)

      data_plot <- 
        dplyr::select(self$data, criteria | starts_with(xvar) ) |> 
                    tidyr::pivot_longer(cols = -xvar, names_to = "criterion") |>
                    dplyr::rename(xvar = 1)

      xlab <- switch(
        xvar,
        "fraction" = expression(paste("|",beta[lambda[1]],"|",{}[1]/max[lambda[1]],"|",beta[lambda[1]],"|",{}[1],sep="")),
        "df" = "Estimated degrees of freedom",
        ifelse(log_scale,expression(log[10](lambda[1])),expression(lambda[1]))
        )
      
      d <- ggplot(data_plot) +
        aes(x=xvar, y=value, colour=criterion, group=criterion) +
        geom_line() + geom_point() + theme_bw() + 
        labs(x = xlab, y = "criterion's value",  title = title)
      
      if (log_scale & (xvar == "lambda")) {
        d <- d + scale_x_log10()
      }
      d
    }
  )
)  
  