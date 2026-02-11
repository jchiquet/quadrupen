context("Consistency of the Lasso solution paths (package 'lars' and 'glmnet')")

tol <- 1e-5

test_that("lasso_quad2lars", {

  require(lars)

  get.lars <- function(x,y,intercept,normalize) {
      lasso.larsen <- lars(x,y,intercept=intercept,normalize=normalize)
      iols <- nrow(lasso.larsen$beta) ## remove last entry corresponding to the OLS estimator
      lambda1 <-  lasso.larsen$lambda ## usde the lars lambda grid
      lasso.quadru <- elastic.net(x,y, intercept=intercept, normalize=normalize,
                                  lambda1=lambda1, lambda2=0, control=list(method="quadra"))
      quad <- list(coef   = as.matrix(lasso.quadru$coefficients),
                   rss    = deviance(lasso.quadru))

      lars <- list(coef   = lasso.larsen$beta[-iols, ],
                   rss    = lasso.larsen$RSS[-iols])

      return(list(quad=quad,lars=lars))
  }

  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)

  ## Run the tests...
  with.intercept <-get.lars(x,y,TRUE,TRUE)
  expect_equal(with.intercept$quad,
              with.intercept$lars, check.attributes = FALSE, tolerance = tol))

  with.intercept.unnormalized <-get.lars(x,y,TRUE,FALSE)
  expect_equal(with.intercept.unnormalized$quad,
              with.intercept.unnormalized$lars, check.attributes = FALSE, tolerance = tol))

  without.intercept <-get.lars(x,y,FALSE,TRUE)
  expect_equal(without.intercept$quad,
              without.intercept$lars, check.attributes = FALSE, tolerance = tol))

  without.intercept.unnormalized <-get.lars(x,y,FALSE,FALSE)
  expect_equal(without.intercept.unnormalized$quad,
              without.intercept.unnormalized$lars, check.attributes = FALSE, tolerance = tol))

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
  with.intercept <-get.lars(x,y,TRUE,TRUE)
  expect_equal(with.intercept$coef.quad,
              with.intercept$coef.lars, check.attributes = FALSE, tolerance = tol))

  with.intercept.unnormalized <-get.lars(x,y,TRUE,FALSE)
  expect_equal(with.intercept.unnormalized$coef.quad,
              with.intercept.unnormalized$coef.lars, check.attributes = FALSE, tolerance = tol))

  without.intercept <-get.lars(x,y,FALSE,TRUE)
  expect_equal(without.intercept$coef.quad,
              without.intercept$coef.lars, check.attributes = FALSE, tolerance = tol))

  without.intercept.unnormalized <-get.lars(x,y,FALSE,FALSE)
  expect_equal(without.intercept.unnormalized$coef.quad,
              without.intercept.unnormalized$coef.lars, check.attributes = FALSE, tolerance = tol))

})

test_that("lasso_quad2glmnet", {

  require(glmnet)

  ## SECOND CHECK: compare to glmnet with prescaling of x
  x <- matrix(rnorm(100*50),100,50)
  y <- rnorm(100)
  y <- y-mean(y)
  n <- nrow(x)
  p <- ncol(x)

  ## If thresh is set to the default, the test won't pass!!!
  ## This is beacause coordinate descent is fast yet not extremely accurate
  lasso.glmn <- glmnet(x,y, lambda.min.ratio=1e-2, thresh=1e-20)
  lasso.quad <- elastic.net(x,y, lambda1=lasso.glmn$lambda*sqrt(n), lambda2=0)

  quad <- list(coef   = as.matrix(lasso.quad$coefficients),
               mu     = lasso.quad$interceptTerm,
               fitted = as.matrix(lasso.quad$fitted))

  glmn <- list(coef   = as.matrix(t(lasso.glmn$beta)),
               mu     = lasso.glmn$a0,
               fitted = predict(lasso.glmn,x))


  expect_equal(quad, glmn, check.attributes = FALSE, tolerance = tol)

})
