context("Consistency of penscale (vs theoretical and 'glmnet')")

tol <- 1e-5

test_that("weighted_ quad2theo", {

  ## SIMPLE CHECK: used identity for design matrix
  n <- 100
  p <- 100
  x <- diag(rep(1,n))
  y <- rnorm(100)

  ## no penscale...
  lasso.quad <- elastic.net(x, y, intercept=FALSE, lambda2=0)
  theo.path  <- t(sapply(lasso.quad$major_tuning, function(lambda) y*pmax(0,1-lambda/abs(y))))
  expect_equal(as.matrix(lasso.quad$coefficients), theo.path, check.attributes = FALSE, tolerance = tol)

  ## glmnet is not equal to what is expected... probably due to the intercept treatment
  lasso.glmn <- glmnet::glmnet(x,y, intercept=FALSE,lambda.min.ratio=1e-2, thresh=1e-20)
  theo.path  <- t(sapply(lasso.glmn$lambda*sqrt(n), function(lambda) y*pmax(0,1-lambda/abs(y))))
  expect_equal(as.matrix(t(lasso.glmn$beta)), theo.path, check.attributes = FALSE, tolerance = 1e-2)

  ## random penscale...
  w <- 1/runif(p,0.5,1)
  w <- w/sum(w)*p ## to fit glmnet rescaling
  lasso.quad <- elastic.net(x, y,intercept=FALSE, penscale = w, lambda2=0)
  theo.path  <- t(sapply(lasso.quad$major_tuning, function(lambda) y*pmax(0,1-lambda*w/abs(y))))
  expect_equal(as.matrix(lasso.quad$coefficients), theo.path, tolerance = 1e-6, check.attributes = FALSE)

  ## glmnet with intercept fit with quadrupen and the theory
  w <- 1/runif(p,0.5,1)
  w <- w/sum(w)*p ## to fit glmnet rescaling
  lasso.glmn <- glmnet::glmnet(x,y, penalty.factor= w,lambda.min.ratio=1e-2, thresh=1e-20)
  lasso.quad <- elastic.net(x,y, lambda1=lasso.glmn$lambda*sqrt(n), penscale = w, lambda2=0)
  
  expect_equal(as.matrix(t(lasso.glmn$beta)), as.matrix(lasso.quad$coefficients), tolerance = tol, check.attributes = FALSE)
  
  ## Check the intercept term also
  expect_equal(lasso.glmn$a0, lasso.quad$intercept, tolerance = tol, check.attributes = FALSE)

})

