## The data is assumed to be centred.


r <- function(A, B, C = NULL, tol = 1e-07) {
  if (!is.null(C)) {
    z <- .Call(stats:::C_Cdqrls, C, cbind(A,B), tol, FALSE)
    A <- as.vector(z$residuals[,1])
    B <- as.vector(z$residuals[,2])
  }
  sum(A * B) / norm(A, "2") / norm(B, "2")
}



r2 <- function(A, B, C = NULL, tol = 1e-07) {
  B <- as.matrix(B)
  if (!is.null(C)) {
    C <- as.matrix(C)
    z <- .Call(stats:::C_Cdqrls, C, cbind(A,B), tol, FALSE)
    A <- as.vector(z$residuals[,1])
    B <- as.matrix(z$residuals[,2:(dim(B)[2]+1)])
  }
  z <- .Call(stats:::C_Cdqrls, B, A, tol, FALSE)
  1 - norm(as.vector(z$residuals), "2")^2 / norm(A, "2")^2
}
