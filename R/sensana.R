
## Creating sensana object
#' @export
sensana <- function(y, d, indep_x, dep_x, x = NULL, z = NULL,
                    quantile = c("normal", "t"), alpha = 0.05) {
  
  quantile <- match.arg(quantile)
  check_sensana(y, d, indep_x, dep_x, x, z, quantile, alpha)
  
  cl <- match.call()
  
  ## Centring and disentangling covariates
  yc <- as.vector(scale(y, scale = FALSE))
  dc <- as.vector(scale(d, scale = FALSE))
  zc <- if(is.null(z)) NULL else as.vector(scale(z, scale = FALSE))
  
  if (is.null(x) || is.null(dep_x)) {
    xtc <- NULL
  } else {
    xtc <- as.matrix(scale(x[, dep_x], scale = FALSE))
    if (length(dep_x) == 1) {colnames(xtc) <- dep_x}
  }
  if (is.null(x) || is.null(indep_x)) {
    xpc <- NULL
  } else {
    xpc <- as.matrix(scale(x[, indep_x], scale = FALSE))
    if (length(indep_x) == 1) {colnames(xpc) <- indep_x}
  }
  
  xc <- cbind(xpc, xtc)
  xzc <- cbind(xc, zc)
  
  
  ## OLS estimate and confidence interval
  if (is.null(xzc)) {
    y_xz <- yc
    d_xz <- dc
  } else {
    decomp_xz <- resid_cpp(xzc, cbind(yc, dc))
    ## decomp_xz <- resid_cpp(cbind(yc, dc), xzc)
    y_xz <- as.vector(decomp_xz[,1])
    d_xz <- as.vector(decomp_xz[,2])
    ## decomp_xz <- .Call(stats:::C_Cdqrls, xzc, cbind(yc, dc), 1e-7, FALSE)
    ## y_xz <- as.vector(decomp_xz$residuals[,1])
    ## d_xz <- as.vector(decomp_xz$residuals[,2])
  }
  y_xzd <- as.vector(resid_cpp(cbind(xzc,dc), as.matrix(yc)))
  ##y_xzd <- as.vector(resid_cpp(as.matrix(yc), cbind(xzc,dc)))
  ## y_xzd <- .Call(stats:::C_Cdqrls, cbind(xzc,dc), yc, 1e-7, FALSE)$residuals
  df <- length(yc) - dim(xzc)[2] - 2 ## intercept + treatment/instrument
  q <- if (quantile == "normal") qnorm(1 - alpha/2) else qt(1 - alpha/2, df)
  
  beta_ols <- sum(y_xz * d_xz) / sum(d_xz^2)
  se_ols <- norm(y_xzd, "2") / norm(d_xz, "2") / sqrt(df)
  confint_ols <- beta_ols + c(-1,1) * q * se_ols
  
  
  ## TSLS estimate and confidence interval
  if (is.null(zc)) {
    beta_tsls <- NULL
    se_tsls <- NULL
    confint_tsls <- NULL
  } else {
    if (is.null(xc)) {
      y_x <- yc
      d_x <- dc
      z_x <- zc
    } else {
      decomp_x <- resid_cpp(xc, cbind(yc, dc, zc))
      ## decomp_x <- resid_cpp(cbind(yc, dc, zc), xc)
      y_x <- as.vector(decomp_x[,1])
      d_x <- as.vector(decomp_x[,2])
      z_x <- as.vector(decomp_x[,3])
      ## decomp_x <- .Call(stats:::C_Cdqrls, xc, cbind(yc, dc, zc), 1e-7, FALSE)
      ## y_x <- as.vector(decomp_x$residuals[,1])
      ## d_x <- as.vector(decomp_x$residuals[,2])
      ## z_x <- as.vector(decomp_x$residuals[,3])
    }
    
    beta_tsls <- sum(y_x * z_x) / sum(d_x * z_x)
    se_tsls <- norm(y_x - beta_tsls*d_x, "2") * norm(z_x, "2") / sum(d_x * z_x) / sqrt(df)
    confint_tsls <- beta_tsls + c(-1,1) * q * se_tsls
  }
    
  ## Initialization of bounds data frame
  bounds <- data.frame(matrix(nrow = 0, ncol = 7))
  colnames(bounds) <- c("arrow", "kind", "lb", "ub", "b", "I", "J")
  
  sa <- list(call = cl,
             y = yc,
             d = dc,
             xt = xtc,
             xp = xpc,
             z = zc,
             beta_ols = beta_ols,
             beta_tsls = beta_tsls,
             confint_ols = confint_ols,
             confint_tsls = confint_tsls,
             se_ols = se_ols,
             se_tsls = se_tsls,
             bounds = bounds,
             alpha = alpha,
             rbc = 0 ## running bound counter
  )
  class(sa) <- "sensana"
  sa
}



#' @export
print.sensana <- function(sa, digits = max(3L, getOption("digits") - 3L)) {
  list2env(sa, environment())
  cat("Sensitivity Analysis:\n\n")
  cat("Dependent Covariates: ", colnames(xt), "\n")
  cat("Independent Covariates: ", colnames(xp), "\n\n")
  
  cat("Estimators:\n")
  cat("OLS\t", round(beta_ols, digits), "\n")
  cat("TSLS\t", round(beta_tsls, digits), "\n\n")
  
  cat(round((1-alpha)*100, digits), "% Confidence Intervals:\n", sep = "")
  cat("OLS\t[", round(confint_ols[1], digits), ",",
      round(confint_ols[2], digits), "]\n")
  cat("TSLS\t[", round(confint_tsls[1], digits), ",",
      round(confint_tsls[2], digits), "]\n\n")
  
  cat("Specified Bounds:\n")
  print(bounds, digits = digits)
  
  cat("\n")
  invisible(sa)
}
