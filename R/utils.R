

## The data is assumed to be centred.
r <- function(A, B, C = NULL) {
  if (!is.null(C)) {
    z <- .Call(stats:::C_Cdqrls, C, cbind(A,B), 1e-07, FALSE)
    A <- as.vector(z$residuals[,1])
    B <- as.vector(z$residuals[,2])
  }
  sum(A * B) / norm(A, "2") / norm(B, "2")
}



r2 <- function(A, B, C = NULL) {
  B <- as.matrix(B)
  if (!is.null(C)) {
    C <- as.matrix(C)
    z <- .Call(stats:::C_Cdqrls, C, cbind(A,B), 1e-07, FALSE)
    A <- as.vector(z$residuals[,1])
    B <- as.matrix(z$residuals[,2:(dim(B)[2]+1)])
  }
  z <- .Call(stats:::C_Cdqrls, B, A, 1e-07, FALSE)
  1 - norm(as.vector(z$residuals), "2")^2 / norm(A, "2")^2
}


finv <- function(f) {
  f / sqrt(1 + f^2)
}


f <- function(r) {
  r / sqrt(1 - r^2)
}

