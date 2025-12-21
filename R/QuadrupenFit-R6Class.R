#' Class "QuadrupenFit"
#'
#' Class of object returned by any fitting function of the
#' \pkg{quadrupen} package (\code{elastic.net} or
#' \code{bounded.reg}).
#' 
#' This class comes with the usual [predict()], [fitted()], [coef()],
#' [residuals()], [show()], [print()] and [deviance()] S3 methods.
#'
#' Specific R6 methods are available for model extraction [QuadrupenFit$get_model()], 
#' cross validation [QuadrupenFit$cross_validate()], stability selection [QuadrupenFit$stability_path()], criteria derivation [QuadrupenFit$criteria()] 
#' and plotting [QuadrupenFit$plot()]. They come with equivalent S3 methods : [cross_validate()], 
#' [stability()] and [plot()].
#'
#' @seealso See also \code{\link{plot.quadrupen}}.
#'
#' @importFrom stats fitted predict residuals deviance
#' @export
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
    infocrit    = NA        ,
    crossval    = NA        ,
    stabsel     = NA        ,
    monitoring  = list()
  ),
  ## ____________________________________________________
  ## 
  ## ACTIVE BINDINGS
  ## ____________________________________________________
  active = list(
    #' @field ncoef number of coefficient (without intercept)
    ncoef = function(value) {private$data$d},
    #' @field nsample sample size
    nsample = function(value) {private$data$n},
    #' @field dataModel an object with class [`DataModel`] storing the data
    dataModel = function(value) {private$data},
    #' @field has_intercept boolean indicating wether an intercept is included in the model
    has_intercept = function(value) {private$intercept},
    #' @field major_tuning vector of "leading" tuning parameters (either l1, linf or l2)
    major_tuning = function(value) {
      if (missing(value))
        return(private$tuning[[1]])
      else private$tuning[[1]] <- value
    },
    #' @field minor_tuning vector of "minor" tuning parameters (either l1 or l2)
    minor_tuning = function(value) {
      if (missing(value))
        return(private$tuning[[2]])
      else private$tuning[[2]] <- value
    },
    #' @field optim_monitoring list monitoring the optimization
    optim_monitoring = function(value) {
      if (!is.null(private$monitoring$convergence))
        private$monitoring$convergence <- 
          sapply(private$monitoring$convergence, status_to_message)
      private$monitoring
    },
    #' @field optim_config list with low level options used for optimization.
    optim_config = function(value) {private$control},
    #' @field fitted Matrix of fitted values, each column corresponding to a value of \code{lambda1}.
    fitted = function(value) {
      Xs <- Matrix::colScale(private$data$X, private$data$norm_X)
      res <- sweep(tcrossprod(Xs, private$beta),2L,-private$mu,check.margin=FALSE)
      res
    },
    #' @field coefficients Matrix (class `"dgCMatrix"`) of
    #' coefficients with respect to the original input. The number of
    #' rows corresponds the length of \code{lambda1}.
    coefficients         = function(value) {
      dimnames(private$beta) <- 
        list(round(c(private$tuning[[1]]),3), colnames(private$data$X))
      private$beta
    },
    #' @field interceptTerm A vector containing the successive values of the 
    #' (unpenalized) intercept.
    #' Equals to zero if \code{intercept} has been set to `FALSE`.
    interceptTerm        = function(value) {private$mu},
    #' @field residuals Matrix of residuals, each column corresponding to a value of `lambda1`.
    residuals            = function(value) {apply(self$fitted, 2, function(y_hat) private$data$y - y_hat)},
    #' @field deviance the model deviance
    deviance             = function(value) {colSums(self$residuals^2)},
    #' @field degrees_freedom Estimated degree of freedoms for the successive `lambda1`.
    degrees_freedom      = function(value) {private$df + ifelse(self$has_intercept, 1L, 0L)},
    #' @field r_squared vector giving the coefficient of determination as a function of lambda1.
    r_squared            = function(value) {1 - colSums(self$residuals^2) / private$data$rss},
    #' @field information_criteria object with class [`Criteria`] storing various information criteria 
    #' (AIC, BIC, GCV, etc) for the current fit.
    information_criteria = function(value) {private$infocrit},
    #' @field cross_validation object with class [`CrossValidation`] storing output of CV job. 
    #' Only available once method cross_validate has been called.
    cross_validation     = function(value) {private$crossval},
    #' @field stability_path object with class [`StabilityPath`] storing output of stability selection. 
    #' Only available once method $stability has been called.
    stability_path       = function(value) {private$stabsel}
  ),
  
  ## ____________________________________________________
  ## 
  ## PUBLIC MEMBERS
  ## ____________________________________________________
  public  = list(
    
    #' @description Initialize a [`QuadrupenFit`] model
    #' @param data a [`DataModel`] object
    #' @param intercept a logical; should an intercept be included in the mode?
    #' @param regParam a list with two elements, a vector and a scalar, for the regularization
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
    #' @description User friendly print method
    show = function() {
      cat("Linear regression with", self$penalty, "penalizer.\n")
      # cat("Linear regression with", x@penalty, "penalizer,", self$rescaling, "rescaling applied to the coefficients.\n")
      if (self$has_intercept) {
        cat("- number of coefficients:", self$ncoef,"+ intercept\n")
      } else {
        cat("- number of coefficients:", self$ncoef,"(no intercept)\n")
      }
      
      cat("- regularization parameter ", names(self$major_tuning), ": ",
          length(self$major_penalty), " points from ",
          format(max(self$major_tuning), digits = 3)," to ",
          format(min(self$major_tuning), digits = 3),"\n", sep="")
      if (!is.null(self$minor_tuning))
        cat("- penalty parameter ",names(self$minor_tuning),": ", self$minor_tuning, "\n", sep="")
      invisible(self)
    },
    #' @description User friendly print method
    print = function() { self$show() },
    #' @description function performing the optimization
    #' @param control list contrlling the optimization process
    fit = function(control) {
      ## ======================================================
      ## C++ CALL OPTIMIZER
      ## 
      if (control$timer) {cpp.start <- proc.time()}
      out <- private$optimizer(private$data, private$tuning, control)
      timer <- ifelse(control$timer, (proc.time() - cpp.start)[3], NA) 
      ## END OF CALL
      ## ======================================================
      private$tuning      <- out$tuning_param
      private$mu          <- drop(out$mu)
      private$beta        <- out$beta
      private$df          <- drop(out$df)
      private$monitoring  <- out$monitoring
      private$monitoring$timer <- timer
      private$control     <- control
    },
    #' @description Model extraction
    #' @param selection either a character (model selection criteria) of a scalar (lambda value)
    #' @param type character for the desired output
    #' @return either a vector of coefficients, a scalar or the model index
    get_model = function(
          selection,
          type = c("coefficients", "penalty", "index")) {

      lambda <- private$tuning[[1]]
      if (is.character(selection)) {
        stopifnot("must be a character in" = selection %in% c("AIC", "BIC", "mBIC", "eBIC", "GCV", "CV_min", "CV_1se"))
        if (selection %in% c("CV_min", "CV_1se") & is.na(private$crossval)) 
          stop("Cross-validation has not yet been performed")
        
        index <- 
          switch(selection,
                 "AIC"  = which.min(private$infocrit$data$AIC),
                 "BIC"  = which.min(private$infocrit$data$BIC),
                 "mBIC" = which.min(private$infocrit$data$mBIC),
                 "eBIC" = which.min(private$infocrit$data$eBIC),
                 "GCV"  = which.min(private$infocrit$data$GCV),
                 "CV_min" = min(match(private$cv_job$lambda1_min, lambda), length(lambda)),
                 "CV_1se" = min(match(private$cv_job$lambda1_1se, lambda), length(lambda)),
          )
      } else {
        index <- match(selection, lambda)
        if (is.na(index)) { ## No exact match
          index <- which.min(abs(selection - lambda)) ## closest model (in terms of parameter value)
          warning(paste("No such a model in the collection. Acceptable parameter values can be found via $major_tuning",
                        paste0("  Returning model with closest value. Requested: ", selection, ", returned: ", lambda[index]),
                        sep = "\n"))
        }
      }
      res <- switch(
        match.arg(type), 
        "index"        = index,
        "penalty"      = lambda[index],
        setNames(c(private$mu[index], private$beta[index, ]), c("intercept", colnames(private$beta)))
        )
      res
    },
    #' @description Predict response for new sample based on the current model
    #' @param newx matrix of new values for the regressor with which to predict. If omitted, the fitted values are used.
    #' @param selection either a character (model selection criteria) of a scalar (lambda value)
    #' @return a vector of predicted value
    predict = function(newx = NULL, selection = NULL) {
      
      if (is.null(selection)) {
        index <- 1:length(private$tuning[[1]])
      } else {
        index <- self$get_model(selection, type = "index")
      }

      if (is.null(newx)) {
        res <- self$fitted[ , index, drop = FALSE]
      } else {
        res <- sweep(newx %*% t(private$beta[index, , drop = FALSE]), 2L, -private$mu[s])
      }
      res
    },
    #' Debias model coefficients
    #' 
    #' @description Apply various debiasing schemes to correct effect of shrinkage 
    #' on the estimation of the coefficients
    #' @param type a character, either "rescaled", "relaxed" or "original". See details
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
      stop("not implemented yet")
    },
    #' Cross-validation for Quadrupen object
    #' 
    #' @description Function that computes K-fold cross-validated error of a
    #' \code{quadrupen} fit, possibly on a grid of `lambda1`, `lambda2`.
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
    #' @return an object with class [`CrossValidation`] is sent back and stored as a 
    #' field of the original [`QuadrupenFit`] object.
    #' 
    #' @seealso [cross_validate()]
    cross_validate = 
      function(
          K       = 10,
          folds   = split(sample(1:self$nsample), rep(1:K, length=self$nsample)),
          lambda2 = self$minor_tuning, verbose = TRUE, cores = detectCores() - 2) {

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
                   mc.cores = cores,
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
    #' Stability selection for Quadrupen object
    #' 
    #' @description Compute the stability path of a (possibly randomized) fitting
    #' procedure as introduced by Meinshausen and Buhlmann (2010).
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
    #' @seealso [stability()]
    stability = function(
          n_subsamples   = 50,
          subsample_size = floor(self$nsample/2),
          subsamples     = replicate(n_subsamples, sample(1:self$nsample, subsample_size), simplify=FALSE),
          weakness       = 1,
          verbose        = TRUE,
          cores       = detectCores() - 2) {
      
      ## =============================================================
      ## INITIALIZATION & PARAMETERS RECOVERY
      if (Sys.info()[['sysname']] == "Windows") {
        warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
        cores <- 1
      }
      
      control  <- private$control; 
      control$verbose <- 0
      control$max.feat  <- self$ncoef
      nlambda1 <- length(private$tuning[[1]])

      ## Prepare blocs of sub samples to run jobs parallely
      blocs <- suppressWarnings(split(1:n_subsamples, 1:cores))
      
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
      prob_bloc <- mclapply(blocs, bloc_stability, mc.cores=cores)
      
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
    #' Penalized criteria based on estimation of degrees of freedom
    #'
    #' @description Produce a plot or send back the values of some penalized criteria
    #' accompanied with the vector(s) of parameters selected
    #' accordingly. The default behavior plots the BIC and the AIC (with
    #' respective factor \eqn{\log(n)}{log(n)} and \eqn{2}{2}) yet the user can specify any
    #' penalty.
    #'
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
    #' @seealso [criteria()]
    criteria = function(penalty=
                          setNames(c(2, log(self$nsample), log(self$ncoef), log(self$nsample) + 2*log(self$ncoef)),
                                   c("AIC","BIC", "mBIC", "eBIC")), sigma=NULL) {
      
      betas <- private$beta
      lambda <- self$major_tuning
      
      n <- self$nsample
      p <- self$ncoef
      
      ## Compute all the penalized criteria
      if (is.null(sigma)) {
        crit <- sapply(penalty, function(pen) n*log(self$deviance/n) + pen * self$degrees_freedom)
      } else {
        crit <- sapply(penalty, function(pen) self$deviance/sigma^2 + pen * self$degrees_freedom)
      }
      crit <- as.data.frame(crit)
      ## Compute generalized cross-validation
      crit$GCV <- self$deviance/(n*(1 - self$degrees_freedom/n)^2)
      
      ## Put together all relevant information about those criteria
      private$infocrit <- 
        Criteria$new(
          value = data.frame(
            crit, 
            df        = self$degrees_freedom, 
            lambda    = lambda, 
            fraction  = apply(abs(betas),1,sum)/max(apply(abs(betas),1,sum)), 
            row.names = 1:nrow(crit)
          )
        )
      invisible(private$infocrit)
    },
    #' Plot method
    #' 
    #' @description Produce a plot of the solution path of a \code{quadrupen} fit.
    #'
    #' @param xvar variable to plot on the X-axis: either \code{"lambda"}
    #' (\eqn{\lambda_1}{lambda1} penalty level or
    #' \eqn{\lambda_2}{lambda2} for ridge regression) or
    #' \code{"fraction"} (\eqn{\ell_1}{l1}-norm
    #' of the coefficients). Default is set to \code{"lambda"}.
    #' @param main the main title. Default is set to the model name followed
    #' by what is on the Y-axis.
    #' @param log.scale logical; indicates if a log-scale should be used
    #' when `xvar="lambda"`. Default is `TRUE`.
    #' @param standardize logical; standardize the coefficients before
    #' plotting (with the norm of the predictor). Default is `TRUE`.
    #' @param labels vector indicating the names associated to the plotted
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
    plot_path = function(xvar = "lambda",
                    main = paste(self$penalty," path", sep=""),
                    log.scale = TRUE, standardize=TRUE, labels = NULL) {
      
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
      print(d)
      
      invisible(d)
    }
  )
)
