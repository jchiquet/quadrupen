library(inline)
library(RcppArmadillo)

code <- '
  arma::mat A  = Rcpp::as<arma::mat>(AA) ;
  arma::vec b  = Rcpp::as<arma::vec>(B) ;
  arma::vec x0 = Rcpp::as<arma::vec>(X) ;

  arma::vec r = b - A * x0;
  arma::vec x = x0;
  arma::vec p = r ;
  double rs_old = dot(r,r) ;
  double rs_new = 1 ;
  int i = 0;
  double alpha ;

  double tol = 1e-2;
  arma::mat Ap ;

  while (sqrt(rs_new) > tol & i < 1e1) {
    Ap = A * p;
    alpha = rs_old/dot(p,Ap) ;
    x += alpha * p ;
    r -= alpha * Ap ;
    rs_new = dot(r,r);
    p = r + rs_new/rs_old*p;
    rs_old = rs_new;
    i++;
  }

   Rcpp::Rcout << "\\n ended after " << i << " iterates." << std::endl;

  return Rcpp::wrap(x) ;
'

cg <- cxxfunction(signature(AA="matrix",B="numeric",X="numeric"),
                  plugin="RcppArmadillo",
                  body=code)

X <- matrix(rnorm(2000*1000),2000,1000)
A <- t(X) %*% X
b <- rnorm(1000)
x0 <- rep(0,1000)

print(system.time(x.reg <- solve(A,b)))
print(system.time(x.cg  <- cg(A,b,x0)))

plot(x.reg, col="red")
points(x.cg, col="blue")
