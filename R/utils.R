## ======================================================
## GRAPH TO MATRIX MANIPULATION

chain_graph <- function(p) {
  G <- list(conn = list(), weight = list())
  G$conn[1:(p - 1)] <- as.integer(2:p) 
  G$conn[p] <- as.integer(1)
  G$weight[1:p] <- rep(1.0, p)
  class(G) <- "FusionGraph"
  G
}

conn_from_adj <- function(adjacency_matrix) {
  G <- list(conn = list(), weight = list())
  pairs <- which(upper.tri(adjacency_matrix), arr.ind = TRUE)
  for (idx in 1:nrow(pairs)) {
    i = unname(pairs[idx,][1])
    j = unname(pairs[idx,][2])
    if(adjacency_matrix[i, j] == 1) {
      G$conn[[i]] = j
      G$conn[[j]] = i
      G$weight[[i]] = G$weight[[j]] <- as.numeric(1)
    }
  }
  class(G) <- "FusionGraph"
  G
}
## ======================================================
## GENERATE A GRID OF PENALTY IF NONE HAS BEEN PROVIDED
get.lambda1.l1 <- function(xty,nlambda1,minratio) {
  lmax <- max(abs(xty))
  return(10^seq(log10(lmax), log10(minratio*lmax), len=nlambda1))
}

get.lambda1.li <- function(xty,nlambda1,minratio) {
  lmax <- sum(abs(xty))
  return(10^seq(log10(lmax), log10(minratio*lmax), len=nlambda1))
}

ctrl_default <- function(d)
  list(verbose     = 0, # default control options
       timer       = FALSE,
       maxiter     = max(500, d),
       method      = "quadra",
       threshold   = 1e-7,
       monitor     = 0,
       bulletproof = TRUE,
       usechol     = TRUE
  )

ctrl_fused_default <- function(d)
  list(verbose       = 0, # default control options
       adjust        = FALSE,
       timer         = FALSE,
       maxiterout    = 100,
       maxiterin     = 10000,
       maxactivation = 10, 
       accuracy      = 1e-6,
       monitor       = 0,
       fusioncheck   = "all", ## c("all","active","none", "naive", "huber")
       usechol       = TRUE
  )

## ====================================================================
## STANDARDIZE THE PREDICTOR (NEEDED FOR CROSS-VALIDATION PURPOSES)
standardize <- function(x,y,intercept,normalize,penscale,zero=.Machine$double.eps) {

  n <- length(y)
  p <- ncol(x)
  ## ============================================
  ## INTERCEPT AND NORMALIZATION TREATMENT
  if (intercept) {
    xbar <- colMeans(x)
    ybar <- mean(y)
  } else {
    xbar <- rep(0,p)
    ybar <- 0
  }

  ## ============================================
  ## NORMALIZATION
  if (normalize) {
    normx <- sqrt(drop(colSums(x^2)- n*xbar^2))
    if (any(normx < zero)) {
      warning("A predictor has no signal: you should remove it.")
      normx[abs(normx) < zero] <- 1 ## dirty way to handle 0/0
    }
    ## normalizing the predictors...
    x <- sweep(x, 2L, normx, "/", check.margin = FALSE)
    ## xbar is scaled to handle internaly the centering of X for
    ## sparsity purpose
    xbar <- xbar/normx
  } else {
    normx <- rep(1,p)
  }
  normy <- sqrt(sum(y^2))

  ## and now normalize predictors according to penscale value
  if (any(penscale != 1)) {
    x <- sweep(x, 2L, penscale, "/", check.margin=FALSE)
    xbar <- xbar/penscale
  }
  ## Computing marginal correlation
  if (intercept) {
    xty   <- drop(crossprod(y-ybar,sweep(x,2L,xbar)))
  } else {
    xty   <- drop(crossprod(y,x))
  }

  return(list(xbar=xbar, ybar=ybar, normx=normx, normy=normy, xty=xty, x=x))
}

status_to_message <- function(status) {
  message <- switch(as.character(status),
                    "0"  = "converged",
                    "1"  = "max # of iterate reached",
                    "2"  = "max # of feature reached",
                    "3"  = "system has become singular",
                    "Return status not recognized"
  )
  message
}

