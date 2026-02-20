context("Consistency of the Elastic-net solution path (package 'elasticnet')")

require(elasticnet)

get.enet <- function(x,y,intercept,normalize=TRUE,method="quadra") {
  lambda2 <- runif(1,0,10)
  enet.larsen <- enet(x,y,lambda=lambda2,intercept=intercept,normalize=normalize)
  iols <- length(enet.larsen$penalty)
  lambda1 <- enet.larsen$penalty[-iols]/2
  
  enet.quadru <- elastic.net(x,y,intercept=intercept,normalize=normalize,
                             lambda1=lambda1, lambda2=lambda2,
                             control = list(method=method))
  
  quad <- as.matrix(enet.quadru$coefficients)
  
  enet <- predict(enet.larsen, type="coefficients", naive=TRUE)$coefficients[-iols,]
  
  return(list(quad = quad, enet = enet))
}

test_that("Elastic-net is correct w.r.t a reference solution", {

  tol <- 1e-2
  
  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)

  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat(" #seed=",seed)
  set.seed(seed)

  beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
  n <- 100
  p <- length(beta)

  mu <- 3 # intercept
  sigma <- 30 # huge noise
  Sigma <- matrix(0.95,p,p) # huge correlation
  diag(Sigma) <- 1

  x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
  y <- 10 + x %*% beta + rnorm(n,0,10)

  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE)
  expect_equal(with.intercept$quad,
              with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  without.intercept <- get.enet(x,y,intercept=FALSE)
  expect_equal(without.intercept$quad,
              without.intercept$enet, check.attributes = FALSE, tolerance = tol)

})


test_that("Elastic-net is correct w.r.t a reference solution - FISTA", {

  tol <- 1e-2
  
  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)
  
  ## Run the tests...
  with.intercept <- get.enet(x,y,intercept=TRUE, method="fista")
  expect_equal(with.intercept$quad,
               with.intercept$enet, check.attributes = FALSE, tolerance = tol)
  
  with.intercept <- get.enet(x,y,intercept=TRUE,normalize=FALSE, method="fista")
  expect_equal(with.intercept$quad,
               with.intercept$enet, check.attributes = FALSE, tolerance = tol)

  ## Run the tests...
  without.intercept <- get.enet(x,y,intercept=FALSE, method="fista")
  expect_equal(without.intercept$quad,
               without.intercept$enet, check.attributes = FALSE, tolerance = tol)
  
  without.intercept <- get.enet(x,y,intercept=FALSE,normalize=FALSE, method="fista")
  expect_equal(without.intercept$quad,
               without.intercept$enet, check.attributes = FALSE, tolerance = tol)
  

})
