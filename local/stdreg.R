varcor<-function(x) {
  l<-lapply(VarCorr(x),function(a) {
    attr(a,"stddev")<-NULL
    attr(a,"correlation")<-NULL
    as.matrix(a)})
  l$sigma2<-sigma(x)^2
  l
}

SS<-function(x) {
  sum((x-mean(x))^2)
}


### This produced standardized y with given betas and X
### X must be standardized
### the code is adapted from

### @MISC {508107,
### TITLE = {Simulate regression with specified standardized coefficients},
### AUTHOR = {whuber (https://stats.stackexchange.com/users/919/whuber)},
### HOWPUBLISHED = {Cross Validated},
### NOTE = {URL:https://stats.stackexchange.com/q/508107 (version: 2021-02-05)},
### EPRINT = {https://stats.stackexchange.com/q/508107},
### URL = {https://stats.stackexchange.com/q/508107}
### }

yield_y<-function(X,beta) {

  n<-dim(X)[1]
  .X<-X/sqrt(n-1)

  dual <- function(aX, tol=1e-16) {
      n <- dim(aX)[1]
      p <- dim(aX)[2]
      qr.Q(qr(cbind(aX, rep(1, n), diag(1, n, n)), tol=tol))[, -seq_len(p+1)]
    }
  y.0 <- .X %*% beta
  y2 <- sum(y.0^2)

  if (y2 > 1) stop("No solutions are possible: betas are too large.")
  eps.norm <- sqrt(1 - y2) # This will be the length of any error term
  Phi <- with(svd(.X), dual(u[, which(d > 1e-8), drop=FALSE]))
##############################################
# This is the entire algorithm.              #
# It uses precomputed y.0, Phi, and eps.norm.#
##############################################
eps <- Phi %*% rnorm(dim(Phi)[2])            # Unnormalized errors
y <- y.0 + eps * eps.norm / sqrt(sum(eps^2)) # Normalize and add
#
# Check the generated `y` meets all constraints.
# maybe we do something with it later on
#  ok=(1 - sqrt(sum(y^2)) > tol*1e2)
#

y*sqrt(n-1)
}
