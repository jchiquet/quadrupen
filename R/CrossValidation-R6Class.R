#' Class CrossValidation
#' 
#' Class of object returned by a cross-validation performed through
#' the \code{$cross_validate()} method.
#' 
#' \@param lambda1 vector of \eqn{\lambda_1}{lambda1}
#' (\eqn{\ell_1}{l1} or \eqn{\ell_\infty}{infinity} penalty levels)
#' for which each cross-validation has been performed.
#' @param lambda2 vector (or scalar) of \eqn{\ell_2}{l2}-penalty levels for
#' which each cross-validation has been performed.
#' @param lambda1.min level of \eqn{\lambda_1}{lambda1} that minimizes the
#' error estimated by cross-validation.
#' @param lambda1.1se largest level of \eqn{\lambda_1}{lambda1} such as
#' the cross-validated error is within 1 standard error of the
#' minimum.
#' @param lambda2.min level of \eqn{\lambda_2}{lambda2} that minimizes the
#' error estimated by cross-validation.
#' @param lambda2.1se largest level of \eqn{\lambda_2}{lambda2}
#'the cross-validated error is within 1 standard error of the
#' minimum (only relevant for ridge regression).
#' @param cv.error a data frame containing the mean
#' cross-validated error and its associated standard error for each
#' values of \code{lambda1} and \code{lambda2}.
#' @param folds list of \code{K} vectors indicating the folds
#' used for cross-validation.
#' 
#' @importFrom dplyr filter
#' 
#' @export
CrossValidation <- R6::R6Class(
  classname = "CrossValidation",
  private = list(
    cv_min = NA,
    cv_1se = NA,
    value  = NA
  ),
  active = list(
    data        = function(value) {private$value},
    lambda1     = function(value) unique(self$error$lambda1),
    lambda2     = function(value) unique(self$error$lambda2),
    lambda1_min = function(value) max(private$cv_min$lambda1),
    lambda2_min = function(value) max(private$cv_min$lambda2),
    lambda1_1se = function(value) max(private$cv_1se$lambda1),
    lambda2_1se = function(value) max(private$cv_1se$lambda2)
  ),
  public = list(
    ## model-related fields
    folds   = "list",
    initialize = 
      function(cv_error, folds) {
        private$value <- cv_error 
        self$folds <- folds 
        private$cv_min <- cv_error |> dplyr::filter(mean <= min(mean))
        private$cv_1se <- cv_error |> 
          dplyr::filter(lambda2 == self$lambda2_min) |> 
          dplyr::filter(mean < min(mean + serr) + 1e-5)
    },
    plotCV_1D = function(main = "Cross-validation error", log.scale=TRUE) {
      ## SIMPLE CROSS-VALIDATION GRAPH
      if (length(self$lambda1) > length(self$lambda2)) {
        d <- ggplot(self$error, aes(x=.data$lambda1,y=.data$mean)) + theme_bw()
        xlab <- ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1]))
        lambda <- data.frame(xval=c(self$lambda1_min,self$lambda1_1se),
                             lambda.choice=factor(c("min. MSE","1-se rule")))
        
      } else {
        ### TODO
        ## ridge or not (meaning working on lambda1 or lambda2)
        d <- ggplot(self$error, aes(x=.data$lambda2,y=.data$mean))
        xlab <- ifelse(log.scale,expression(log[10](lambda[2])),expression(lambda[2]))
        lambda <- data.frame(xval=c(self$lambda2_min,self$lambda2_1se),
                             lambda.choice=factor(c("min. MSE","1-se rule")))
      }
      d <- d + xlab(xlab) + ylab("Mean square error") + 
        geom_point(alpha=0.3) + 
        geom_smooth(aes(ymin=.data$mean-.data$serr, ymax=.data$mean+.data$serr, group=as.factor(.data$lambda2),
                        colour=as.factor(.data$lambda2)), data=self$error, alpha=0.2, stat="identity")
      if (log.scale) {
        d <- d + scale_x_log10() + annotation_logticks(sides="b")
      }
      d <- d + ggtitle(main) +
        geom_vline(data=lambda, aes(xintercept=.data$xval, linetype=.data$lambda.choice),
                   alpha=0.2, show.legend = TRUE) + 
          scale_color_discrete(labels = round(self$lambda2, 3)) +
        labs(colour = expression(explored~lambda[2]), linetype = expression(lambda[1]~choice))
      
      d
    },
    plotCV_2D = function(main = "Cross-validation error") {
      d <- ggplot(data=self$error, aes(x=.data$lambda1, y=.data$lambda2, z=.data$mean))
      d <- d + geom_tile(aes(fill=.data$mean)) + stat_contour(linewidth=0.2, binwidth=diff(range(self$error$mean))/10) + ggtitle(main)
      d <- d + scale_x_continuous(trans=log10_trans()) + xlab(expression(log[10](lambda[1])))
      d <- d + scale_y_continuous(trans=log10_trans()) + ylab(expression(log[10](lambda[2])))
      d <- d + annotation_logticks() + theme_bw()
      i_1se <- which(self$error$mean - self$error$serr <= min(self$error$mean))
      d <- d + stat_contour(alpha=0.5, colour="#CCCCCC", linewidth=0.65, breaks=quantile(self$error$mean[i_1se], probs=seq(0,1,len=6)))
      d
    },
    plot = function(log.scale=TRUE, plot=TRUE, main = "Cross-validation error", ...) {
      if (length(self$lambda1) > 1 & length(self$lambda2) > 1) {
        d <- self$plotCV_2D(main)
      } else {
        d <- self$plotCV_1D(main, log.scale)  
      }
      ## DO THE PLOT
      if (plot) print(d)
      invisible(d)
    }
  )
)
