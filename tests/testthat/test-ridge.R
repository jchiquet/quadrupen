context("Consistency of the Lasso solution paths (package 'lars' and 'glmnet')")

tol <- 1e-4
require(MASS)

test_that("lasso_quad2lars", {
  
  
  get_ridge <- function(x,y,intercept) {

    lambda <- 10^seq(2,-2,len = 100)
    if (intercept) {
      ridge_mass <- MASS::lm.ridge(y ~ x + 1, lambda = lambda * nrow(x))
      mass <- list(coef = t(coefficients(ridge_mass)[, -1]),
                   mu   = coefficients(ridge_mass)[, 1])
    } else {
      ridge_mass <- MASS::lm.ridge(y ~ x + 0, lambda = lambda * nrow(x))
      mass <- list(coef = t(coefficients(ridge_mass)),
                   mu   = rep(0,length(lambda)))
      
    }
    ridge_quad <- ridge(x, y, intercept=intercept, normalize=TRUE, lambda=lambda)
    quad <- list(coef = as.matrix(ridge_quad$coefficients),
                 mu   = ridge_quad$intercept)
    
    
    return(list(quad=quad,mass=mass))
  }
  
  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)
  
  ## Run the tests...
  with.intercept <- get_ridge(x,y,TRUE)
  expect_equal(with.intercept$quad,
               with.intercept$mass, check.attributes = FALSE, tolerance = tol)
  
  without.intercept <- get_ridge(x,y,FALSE)
  expect_equal(without.intercept$quad,
               without.intercept$mass, check.attributes = FALSE, tolerance = tol)
  
  ## RANDOM DATA
  seed <- sample(1:10000,1)
  ## cat("\n#seed=",seed)
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
  ## Run the tests...
  with.intercept <- get_ridge(x,y,TRUE)
  expect_equal(with.intercept$quad,
               with.intercept$mass, check.attributes = FALSE, tolerance = tol)
  
  without.intercept <- get_ridge(x,y,FALSE)
  expect_equal(without.intercept$quad,
               without.intercept$mass, check.attributes = FALSE, tolerance = tol)
})
