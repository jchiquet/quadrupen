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
QuadrupenFit <- R6::R6Class(
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
    tuning      = numeric() ,
    intercept   = NA        ,
    control     = list()    ,
    optimizer   = NA        ,
    crossval    = NA        ,
    stabsel     = NA        ,
    monitoring  = list()
  ),
  ## ____________________________________________________
  ## 
  ## ACTIVE BINDINGS
  ## ____________________________________________________
  active = list(
    ncoef = function(value) {private$data$d},
    nsample = function(value) {private$data$n},
    dataModel = function(value) {private$data},
    has_intercept = function(value) {private$intercept},
    #' @field major_penalty vector of "leading" tuning parameters (either l1, linf or l2)
    major_tuning = function(value) {
      if (missing(value))
        return(private$tuning[[1]])
      else private$tuning[[1]] <- value
    },
    #' @field minor_penalty vector of "minor" tuning parameters (either l1 or l2)
    minor_tuning = function(value) {
      if (missing(value))
        return(private$tuning[[2]])
      else private$tuning[[2]] <- value
    },
    optim_monitoring = function(value) {private$monitoring},
    optim_config = function(value) {private$control},
    fitted = function(value) {
      Xs <- Matrix::colScale(private$data$X, private$data$norm_X)
      res <- sweep(tcrossprod(Xs, private$beta),2L,-private$mu,check.margin=FALSE)
      res
    },
    coefficients     = function(value) {private$beta},
    interceptTerm    = function(value) {private$mu},
    residuals        = function(value) {apply(self$fitted, 2, function(y_hat) private$data$y - y_hat)},
    deviance         = function(value) {colSums(self$residuals^2)},
    degrees_freedom  = function(value) {private$df + ifelse(self$has_intercept, 1L, 0L)},
    r_squared        = function(value) {1 - colSums(self$residuals^2) / private$data$rss},
    cv_job           = function(value) {private$crossval},
    stabsel_job      = function(value) {private$stabsel}
  ),
  ## ____________________________________________________
  ## 
  ## PUBLIC MEMBERS
  ## ____________________________________________________
  public  = list(
    initialize = function(data, intercept, regParam) {

      stopifnot("The data object must be an instance of DataModel" = inherits(data, "DataModel"))
      stopifnot("regParam must be a list" = is.list(regParam))
      stopifnot("All regularization parameters must be positive." = all(unlist(regParam) >= 0))
      stopifnot("The first entry of regParam must be sorted in decreasing order." =
                  !is.unsorted(rev(regParam[[1]])))
      if (length(regParam) > 1)
        stopifnot("The second entry of regParam must be a scalar (cannot be a vector)." = 
                  (length(regParam[[2]]) == 1 & inherits(regParam[[2]], "numeric")))

      private$data      <- data
      private$intercept <- intercept
      private$tuning    <- regParam
    },
    show = function() {
      cat("Linear regression with", self$penalty, "penalizer.\n")
      # cat("Linear regression with", x@penalty, "penalizer,", self$rescaling, "rescaling applied to the coefficients.\n")
      if (self$has_intercept) {
        cat("- number of coefficients:", self$ncoef,"+ intercept\n")
      } else {
        cat("- number of coefficients:", self$ncoef,"(no intercept)\n")
      }
      
      cat("- regularization parameter ",names(self$major_tuning), ": ",
          length(self$major_penalty), " points from ",
          format(max(self$major_tuning), digits = 3)," to ",
          format(min(self$major_tuning), digits = 3),"\n", sep="")
      cat("- penalty parameter ",names(self$minor_tuning),": ", self$minor_tuning, "\n", sep="")
      invisible(self)
    },
    #' @description User friendly print method
    print = function() { self$show() },
    predict = function(newx = NULL, ... ) {
      if (is.null(newx)) {
        res <- self$fitted
      } else {
        res <- sweep(newx %*% t(private$beta), 2L, -private$mu)
      }
      res
    },
    
    #' Debiasing for quadrupen object
    #' 
    #' @description Apply various debiasing schemes to correct effect of shrinkage 
    #' on the estimation of the coefficients
    #' 
    #' @param type a character, either "rescaled", "relaxed" or "original. See details
    #' 
    #' @details
    #' the 'rescaled' debaising is as defined in Zou and Hastie (2006): the vector of
    #' parameters is rescaled by a coefficient \code{(1+lambda2)}. 'Original' reset to 
    #' the original scaling
    #' 
    #' @return nothing is return, the beta are internaly rescaled
    #'
    debias = function(type = c("rescaled")) {
      ### TODO
    },
    #' Cross-validation for quadrupen object
    #' 
    #' @description Function that computes K-fold cross-validated error of a
    #' \code{quadrupen} fit, possibly on a grid of
    #' \code{lambda1,lambda2}.
    #'
    #' @param K integer indicating the number of folds. Default is 10.
    #'
    #' @param folds list of \code{K} vectors that describes the folds to
    #' use for the cross-validation. By default, the folds are randomly
    #' sampled with the specified K. The same folds are used for each
    #' values of \code{lambda2}.
    #'
    #' @param lambda2 tunes the \eqn{\ell_2}{l2}-penalty (ridge-like) of
    #' the fit. If none is provided, a vector of values is generated and
    #' a CV is performed on a grid of \code{lambda2} and \code{lambda1},
    #' using the same folds for each \code{lambda2}). Ignored when
    #' \code{penalty} equals \code{"lasso"}. CV is only performed on
    #' \code{lambda2} when the \code{ridge} penalty is used.
    #'
    #' @param verbose logical; indicates if the progression (the current
    #' lambda2) should be displayed. Default is \code{TRUE}.
    #'
    #' @param mc.cores the number of cores to use. The default uses all
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
    #' @return An object of class "CrossValidation" for which a \code{plot} method
    #' is available.
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
    #' ## Use fewer lambda1 values by overwritting the default parameters
    #' ## and cross-validate over the sequences lambda1 and lambda2
    #' cv.grid <- crossval(x,y, lambda2=10^seq(2,-2,len=50), nlambda1=50)
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
    cross_validate = 
      function(
          K       = 10,
          folds   = split(sample(1:self$nsample), rep(1:K, length=self$nsample)),
          lambda2 = self$minor_tuning, verbose = TRUE, mc.cores = detectCores() - 2) {

        ## Some variables and copies useful for CV work
        K <- length(folds)
        control  <- private$control; control$verbose <- 0
        regParam <- private$tuning # copy of the currently used tuning parameters
        nlambda1 <- length(regParam[[1]])
        nlambda2 <- length(lambda2)
        lambda2_vec <- rep(lambda2, each = K)
        fold_id <- rep(1:K, nlambda2)
        
        if (verbose){
          cat("\nCROSS-VALIDATION FOR ", self$penalty," REGULARIZER \n\n")
          cat(K, "-fold CV on a grid of (", 
              nlambda1, ",", nlambda2, ") tuning parameters\n", sep="")
        }

        ## Same data splitting is kept for varying lambda2 values
        CVData <- self$dataModel$splitTrainTest(K, folds)
        for (fold in 1:K) {
          CVData[[fold]]$trainData$standardize(self$dataModel$is_centered, self$dataModel$is_scaled)
          CVData[[fold]]$trainData$getSufficientStat()
        }

        ## CV err for a fixed couple fold/lambda2
        one_fold <- function(fold, lambda2_) {
          if (verbose & (fold == 1)) cat(round(lambda2_, 3),"\t")
          regParam[[2]] <- lambda2_
          out <- private$optimizer(CVData[[fold]]$trainData, regParam, control)
          y_hat <- scale(tcrossprod(CVData[[fold]]$testData$X, out$beta), -out$mu, FALSE)
          err <- sweep(y_hat, 1L, CVData[[fold]]$testData$y)^2
          if (ncol(err) < length(regParam[[1]])) {
            NAs <- length(regParam[[1]]) - ncol(err)
            err <- cbind(err, matrix(NA, nrow(err),NAs))
          }
          err
        }

        err <- do.call(rbind, 
          parallel::mcmapply(FUN = one_fold, fold = fold_id, lambda = lambda2_vec, 
                   mc.cores = mc.cores,
                   mc.preschedule = ifelse(K > 10,TRUE,FALSE),
                   SIMPLIFY = FALSE
          )) |> as.matrix() |> as.data.frame()
        if (verbose) cat("\n")
        
        res <- do.call(rbind, tapply(err, rep(1:nlambda2, each=self$nsample), function(err_) {
          mean <- colMeans(err_, na.rm = TRUE)
          if (any(is.nan(mean))) {
            warning("\nThere have been a lot of early stops along the path: 
                  I keep on running, but you should reconsider the value of 
                  the minimal penalty along the path regarding the n<<p setting.")
          }
          mean[is.nan(mean)] <- NA
          serr <- colSums(sweep(err_, 2L, mean, check.margin = FALSE)^2, na.rm = TRUE)
          serr <- sqrt((serr/(self$dataModel$n - 1))/K)
          data.frame(mean = mean, serr = serr, lambda1 = private$tuning[[1]])
        }))
        res$lambda2 <- rep(lambda2, each = length(private$tuning[[1]]))

        private$crossval <- CrossValidation$new(cv_error = res, folds = folds)
        invisible(private$crossval)
    },
    #' @description Compute the stability path of a (possibly randomized) fitting
    #' procedure as introduced by Meinshausen and Buhlmann (2010).
    #'
    #' @param x matrix of features, possibly sparsely encoded
    #' (experimental). Do NOT include intercept.
    #'
    #' @param y response vector.
    #'
    #' @param n_subsamples integer indicating the number of subsamplings
    #' used to estimate the selection probabilities. Default is 100.
    #'
    #' @param subsample_size integer indicating the size of each subsamples.
    #' Default is \code{floor(n/2)}.
    #'
    #' @param subsamples list with \code{subsamples} entries with vectors
    #' describing the folds to use for the stability procedure. By
    #' default, the folds are randomly sampled with the specified
    #' \code{n_subsamples} and \code{subsample_size} argument.
    #'
    #' @param weakness Coefficient used for randomizing the weights of each features.
    #' Default is \code{0.5}. Set to 1 for no randomization. See
    #' details below.
    #' 
    #' @param verbose logical; indicates if the progression should be
    #' displayed. Default is \code{TRUE}.
    #'
    #' @param mc.cores the number of cores to use. The default uses all
    #' the cores available.
    #'
    #' @return An object of class \code{StabilityPath}.
    #'
    #' @note When \code{randomized = TRUE}, the \code{penscale} argument
    #' that weights the penalty tuned by \eqn{\lambda_1}{lambda1} is
    #' perturbed (divided) for each subsample by a random variable
    #' uniformly distributed on
    #' \if{latex}{\eqn{[\alpha,1]}}\if{html}{[&#945;,1]}\if{text}{\eqn{[alpha,1]}},
    #' where
    #' \if{latex}{\eqn{\alpha}}\if{html}{&#945;}\if{text}{\eqn{alpha}} is
    #' the weakness parameter.
    #'
    #' If the user runs the fitting method with option
    #' \code{'bulletproof'} set to \code{FALSE}, the algorithm may stop
    #' at an early stage of the path. Early stops of the underlying
    #' fitting function are handled internally, in the following way: we
    #' chose to simply skip the results associated with such runs, in
    #' order not to bias the stability selection procedure. If it occurs
    #' too often, a warning is sent to the user, in which case you should
    #' reconsider the grid of \code{lambda1} for stability selection. If
    #' \code{bulletproof} is \code{TRUE} (the default), there is nothing
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
    stability = function(
          n_subsamples   = 50,
          subsample_size = floor(self$nsample/2),
          subsamples     = replicate(n_subsamples, sample(1:self$nsample, subsample_size), simplify=FALSE),
          weakness       = 1,
          verbose        = TRUE,
          mc.cores       = detectCores() - 2) {
      
      ## =============================================================
      ## INITIALIZATION & PARAMETERS RECOVERY
      if (Sys.info()[['sysname']] == "Windows") {
        warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
        mc.cores <- 1
      }
      
      control  <- private$control; 
      control$verbose <- 0
      control$max.feat  <- self$ncoef
      nlambda1 <- length(private$tuning[[1]])

      ## Prepare blocs of sub samples to run jobs parallely
      blocs <- suppressWarnings(split(1:n_subsamples, 1:mc.cores))
      
      if (verbose) {
        cat(paste("\n\nSTABILITY SELECTION ",ifelse(weakness < 1,"with","without")," randomization (weakness = ",weakness,")",sep=""))
        cat(paste("\nFitting procedure: ", self$penalty," with ", nlambda1,"-dimensional grid of lambda1.", sep=""))
        cat("\nRunning",length(blocs),"jobs parallely (1 per core)")
        cat("\nApprox.", length(blocs[[1]]),"subsamplings for each job for a total of", n_subsamples)
      }
      ## get data samples
      SubsampledData <- self$dataModel$splitSubSamples(n_subsamples, subsample_size, subsamples, weakness)

      ## function to run on each core
      bloc_stability <- function(subsets) {
        select <- Matrix(0, nlambda1, self$ncoef)
        subsamples_ok <- 0
        for (s in 1:length(subsets)) {
          active <- private$optimizer(SubsampledData[[subsets[s]]], private$tuning, control)$active
          if (nrow(active) == nlambda1) {
            subsamples_ok <- subsamples_ok + 1
            select <- select + active
          }
        }
        if (subsamples_ok < 0.5*length(subsets)) {
          cat("\nWarning: more than 50% of the run were discarded in that core due to early stops of the fitting procedure. You should consider largest 'min.ratio' or strongest 'lambda2'.")
        }
        return(select/(subsamples_ok*length(blocs)))
      }
      
      ## Now launch the B jobs...
      prob_bloc <- mclapply(blocs, bloc_stability, mc.cores=mc.cores)
      
      ## Construct the probability path
      path <- Matrix(0, nlambda1, self$ncoef)
      for (b in 1:length(prob_bloc)) {
        path <- path + prob_bloc[[b]]
      }

      private$stabsel <- StabilityPath$new(
        probabilities = path          ,
        regParam      = private$tuning,
        subsamples    = subsamples
      )
      invisible(private$stabsel)
    },
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
    #'
    plot = function(xvar = "lambda",
                    main = paste(self$penalty," path", sep=""),
                    log.scale = TRUE, standardize=TRUE, labels = NULL, plot = TRUE, ...) {
      
      if (length(self$major_tuning) == 1) {
        stop("Not available when the leading vector of tuning parameters boild down to a scalar.")
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
        xv <- self$major_tuning
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
    criteria = function(penalty=
        setNames(c(2, log(self$nsample), log(self$ncoef), log(self$nsample) + 2*log(self$ncoef)),
                 c("AIC","BIC", "mBIC", "eBIC")), sigma=NULL,
                         log.scale=TRUE, xvar = "lambda", plot=TRUE) {
      
      betas <- private$beta
      lambda <- self$major_tuning
      
      n <- self$nsample
      p <- self$ncoef
      

      ## compute the penalized criteria
      if (is.null(sigma)) {
        crit <- sapply(penalty, function(pen) n*log(self$deviance/n) + pen * self$degrees_freedom)
      } else {
        crit <- sapply(penalty, function(pen) self$deviance/sigma^2 + pen * self$degrees_freedom)
      }
      crit <- as.data.frame(crit)

      ## Compute generalized cross-validation
      crit$GCV <- self$deviance/(n*(1 - self$degrees_freedom/n)^2)
      
      ## recover the associated vectors of parameters
      beta.min  <- t(betas[apply(crit, 2, which.min), ])
      if (!is.null(dim(beta.min)))
        colnames(beta.min) <- c(names(penalty), "GCV")

      if (inherits(self$cv_job, "CrossValidation")) {
        crit$CV <- self$cv_job$error$mean
        i_min <- min(match(self$cv_job$lambda1_min, lambda),length(lambda))
        i_1se <- min(match(self$cv_job$lambda1_1se, lambda),length(lambda))
        beta_CV <- t(betas[c(i_min, i_1se), ])
        colnames(beta_CV)  <- c("CV_min", "CV_1se")
        beta.min <- cbind(beta.min, beta_CV)
      }
      
      ## put together all relevant information about those criteria
      criterion <- data.frame(crit, df=self$degrees_freedom, lambda=lambda, fraction = apply(abs(betas),1,sum)/max(apply(abs(betas),1,sum)), row.names=1:nrow(crit))

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
          geom_line(aes(x=xvar,y=value)) + geom_point(aes(x=xvar,y=value)) + theme_bw() + 
          labs(x=xlab, y="criterion's value",  title=paste("Information Criteria for a", self$penalty,"fit"))
        
        if (log.scale & (xvar=="lambda")) {
          d <- d + scale_x_log10()
        }
        print(d)
        return(invisible(d))
      } else {
        return(list(criterion=criterion, beta.min=beta.min))
      }
    }
  )
)

## Auxiliary functions to check the given class of an object
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
