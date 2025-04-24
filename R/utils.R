
## The data must be centred.

r <- function(A, B, C = NULL) {
  if (!is.null(C)) {
    res <- r_partial(A, B, C)
  } else {
    res <- r_marginal(A, B)
  }
  res
}


r2 <- function(A, B, C = NULL) {
  B <- as.matrix(B)
  if (!is.null(C)) {
    C <- as.matrix(C)
    ret <- r2_partial(A, B, C)
  } else {
    ret <- r2_marginal(A, B)
  }
  ret
}


finv <- function(f) {
  f / sqrt(1 + f^2)
}


f <- function(r) {
  r / sqrt(1 - r^2)
}

