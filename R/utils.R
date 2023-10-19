## Functions for R- and R^2-values
## It's assumed that all data is centred before.
proj <- function(a, weights){
  ##a %*% solve(t(a) %*% a) %*% t(a)
  if (is.null(weights)) {
    ret <- a %*% MASS::ginv(t(a) %*% a, tol = 0) %*% t(a)
  } else {
    w <- diag(weights)
    ret <- a %*% MASS::ginv(t(a) %*% w %*% a, tol = 0) %*%
      t(a) %*% w
  }
  ret
}

r <- function(y, x, z = NULL, weights = NULL) {
  if (!is.null(z)) {
    z <- as.matrix(z)
    pz <- proj(z, weights)
    qz <- diag(ncol(pz)) - pz
    y <- as.vector(qz %*% y)
    x <- as.vector(qz %*% x)
  }
  if (!is.null(weights)) {
    y <- sqrt(weights) * y
    x <- sqrt(weights) * x
  }
  sum(x * y) / norm(x, "2") / norm(y, "2")
}

r2 <- function(y, x, z = NULL, weights = NULL) {
  x <- as.matrix(x)
  if (!is.null(z)) {
    z <- as.matrix(z)
    pz <- proj(z, weights)
    qz <- diag(ncol(pz)) - pz
    y <- as.vector(qz %*% y)
    x <- qz %*% x
  }
  ypx <- as.vector(proj(x, weights) %*% y)
  if (!is.null(weights)) {
    y <- sqrt(weights) * y
    ypx <- sqrt(weights) * ypx
  }
  norm(ypx, "2")^2 / norm(y, "2")^2
}


compute_r2_se <- function(y, d, x, z, df, weights = NULL) {
  p_xzd <- proj(cbind(x, z, d), weights)
  p_xz <- proj(cbind(x, z), weights)
  q_xzd <- diag(ncol(p_xzd)) - p_xzd
  q_xz <- diag(ncol(p_xz)) - p_xz

  y_xzd <- as.vector(q_xzd %*% y)
  d_xz <- as.vector(q_xz %*% d)
  se <- sd(y_xzd) / sd(d_xz) / sqrt(df)

  y_xz <- as.vector(q_xz %*% y)
  beta_ols <- cov(d_xz, y_xz) / var(d_xz)

  c(beta_ols, se)
}
