context("Testing consistency and timings of warm restart")

test_that("warm_restart", {

  require(quadrupen)

  get.coef <- function(x,y) {
    lambda1 <- .25
    
    enet.ref <- elastic.net(x,y,lambda1=lambda1, control=list(timer=TRUE))

    enet.ref.bot <- elastic.net(x,y,lambda1=lambda1*2)
    enet.ref.up  <- elastic.net(x,y,lambda1=lambda1/2)
    beta0_bt <- as.numeric(enet.ref.bot$coefficients)
    beta0_up <- as.numeric(enet.ref.up$coefficients)
    enet.bot <- elastic.net(x,y,lambda1=lambda1, beta0 = beta0_bt, control=list(timer=TRUE))
    enet.up  <- elastic.net(x,y,lambda1=lambda1, beta0 = beta0_up, control=list(timer=TRUE))

    cat("\n\tTimings with warm-restart along the path")
    cat("\n\t\tfrom stratch: ",enet.ref$optim_monitoring$timer)
    cat("\n\t\tstarting from sparser solution: ",enet.bot$optim_monitoring$timer)
    cat("\n\t\tstarting from more dense solution: ",enet.up$optim_monitoring$timer)
    cat("\n\n")
    
    return(list(
      coef.ref=as.matrix(enet.ref$coefficients),
      coef.bot=as.matrix(enet.bot$coefficients),
      coef.up =as.matrix(enet.up$coefficients)))
  }

  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)

  ## Run the tests...
  cat("\n  * tiny-size problem...")
  out <- get.coef(x,y)
  expect_equal(out$coef.bot, out$coef.ref, check.attributs = FALSE)
  expect_equal(out$coef.up , out$coef.ref, check.attributs = FALSE)

  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat(" #seed=",seed)
  set.seed(seed)

  beta <- rep(rep(c(0,1,0,-1,0), c(25,10,25,10,25)),5)
  n <- 300
  p <- length(beta)

  mu <- 3 # intercept
  sigma <- 30 # huge noise
  Sigma <- matrix(0.95,p,p) # huge correlation
  diag(Sigma) <- 1

  x <- as.matrix(matrix(rnorm(p*n),n,p) %*% chol(Sigma))
  y <- 10 + x %*% beta + rnorm(n,0,10)

  ## Run the tests...
  cat("\n  * small-size problem, with correlation...")
  out <- get.coef(x,y)
  # expect_equal(out$coef.bot, out$coef.ref, check.attributs = FALSE)
  # expect_equal(out$coef.up , out$coef.ref, check.attributs = FALSE)

})

