
#' Sensitivity Analysis for Linear Regression and Instrumental Variable models
#' 
#' `sensana()` does data pre-processing for the sensitivity analysis, computes the
#' OLS and TSLS estimates and confidence intervals and creates a
#' `sensana` object. 
#' 
#' 
#' @details
#' We split the covariates \eqn{X} into 2 groups: \eqn{\dot{X}}
#' and \eqn{\tilde{X}}. The arguments `indep_x` and `dep_x` contain the indices/names
#' of the respective covariates. For an unmeasured confounder \eqn{U}, we assume
#' that \eqn{\dot{X} \perp U \mid \tilde{X}} holds.
#' 
#' This method (and package) can be used even if there are no instrumental variables
#' or covariates. The absence of these is indicated by passing `NULL` as an argument
#' to the respective parameters.
#' 
#' To specify a sensitivity model, use the method [add_bound()]. To conduct inference,
#' use the method [pir()] to estimate the partially identified range and use the
#' method [sensint()] to estimate a sensitivity interval.
#' 
#' More details are provided in the article cited below.
#' 
#' 
#' @param y Numeric vector of outcomes.
#' @param d Numeric vector of treatments/endogeneous variables.
#' @param indep_x Character vector of names of cond. independent covariates. Default is `NULL`.
#' @param dep_x Character vector of names of remaining covariates. Default is `NULL`.
#' @param x Numeric matrix or data.frame of covariates. Default is `NULL`.
#' @param z Numeric vector of instruments. Default is `NULL`.
#' @param quantile Choice of quantile used in confidence interval: `"normal"`
#'  (default) or `"t"`.
#' @param alpha Significance level for the confidence interval. Default is `0.05`.
#' 
#' @returns Object of the class `sensana`. This is a list
#'  containing the following components:
#'  \item{call}{The function call.}
#'  \item{y}{The centered entries of the \code{y} vector.}
#'  \item{d}{The centered entries of the \code{d} vector.}
#'  \item{xt}{The centered matrix of covariates indexed by \code{dep_x}.}
#'  \item{xp}{The centered matrix of covariates indexed by \code{indep_x}.}
#'  \item{z}{The centered entries of the \code{z} vector.}
#'  \item{beta_ols}{OLS estimate.}
#'  \item{beta_tsls}{TSLS estimate.}
#'  \item{confint_ols}{\code{1-alpha} OLS confidence interval.}
#'  \item{confint_tsls}{\code{1-alpha} TSLS confidence interval.}
#'  \item{se_ols}{Standard error of OLS estimator.}
#'  \item{se_tsls}{Standard error of TSLS estimator.}
#'  \item{bounds}{An empty data frame that stores bounds on the sensitivity
#'  parameters. Users can add bounds via [add_bound()].}
#'  \item{alpha}{Specified confidence level.}
#'  \item{rbc}{Internal counter used for default names of bounds.}
#' 
#' @examples
#' x <- rnorm(20)
#' z <- 0.5 * x + rnorm(20)
#' d <- z - 0.5 * x + rnorm(20)
#' y <- d + x + rnorm(20)
#' sa <- sensana(y = y, d = d, x = data.frame(x = x), z = z,
#'               indep_x = "x", dep_x = NULL, alpha = 0.05, quantile = "t")
#'               
#'
#' @references
#' Freidling T, Zhao Q (2024). “Optimization-based Sensitivity Analysis for
#' Unmeasured Confounding using Partial
#' Correlations.” _arXiv preprint arXiv:2301.00040v3_.
#' @export
sensana <- function(y, d, indep_x, dep_x, x = NULL, z = NULL,
                    quantile = c("normal", "t"), alpha = 0.05) {
  
  quantile <- match.arg(quantile)
  check_sensana(y, d, indep_x, dep_x, x, z, quantile, alpha)
  
  cl <- match.call()
  
  ## Centering and disentangling covariates
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
    y_xz <- as.vector(decomp_xz[,1])
    d_xz <- as.vector(decomp_xz[,2])
  }
  y_xzd <- as.vector(resid_cpp(cbind(xzc,dc), as.matrix(yc)))

  n <- length(yc)
  df_ols <- n - qr(cbind(xzc, dc, rep(1,n)))$rank
  
  q_ols <- if (quantile == "normal") stats::qnorm(1 - alpha/2) else stats::qt(1 - alpha/2, df_ols)
  
  beta_ols <- sum(y_xz * d_xz) / sum(d_xz^2)
  se_ols <- norm(y_xzd, "2") / norm(d_xz, "2") / sqrt(df_ols)
  confint_ols <- beta_ols + c(-1,1) * q_ols * se_ols
  
  
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
      df_tsls <- n - 2
    } else {
      decomp_x <- resid_cpp(xc, cbind(yc, dc, zc))
      y_x <- as.vector(decomp_x[,1])
      d_x <- as.vector(decomp_x[,2])
      z_x <- as.vector(decomp_x[,3])
      df_tsls <- n - qr(cbind(xc,rep(1,n)))$rank - 1
    }
    
    q_tsls <- if (quantile == "normal") stats::qnorm(1 - alpha/2) else stats::qt(1 - alpha/2, df_tsls)
    
    beta_tsls <- sum(y_x * z_x) / sum(d_x * z_x)
    se_tsls <- norm(y_x - beta_tsls*d_x, "2") * norm(z_x, "2") / sum(d_x * z_x) / sqrt(df_tsls)
    confint_tsls <- beta_tsls + c(-1,1) * q_tsls * se_tsls
  }
    
  ## Initialization of data frame for bounds
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


#' Print Method for `sensana` objects
#'
#' @param x An object of class \code{"sensana"}.
#' @param digits The number of significant digits to display.
#' @param ... Additional arguments passed to or from other methods (ignored here).
#'
#' @method print sensana
#' @return Invisibly returns \code{x}, the original object.
#' @examples
#' x <- rnorm(20)
#' z <- 0.5 * x + rnorm(20)
#' d <- z - 0.5 * x + rnorm(20)
#' y <- d + x + rnorm(20)
#' sa <- sensana(y = y, d = d, x = data.frame(x = x), z = z,
#'               indep_x = "x", dep_x = NULL, alpha = 0.05, quantile = "t")
#' print(sa)
#' @export
print.sensana <- function(x, digits = max(3L, getOption("digits") - 3L),...) {
  cat("Sensitivity Analysis:\n\n")
  cat("Dependent Covariates: ", colnames(x$xt), "\n")
  cat("Independent Covariates: ", colnames(x$xp), "\n\n")
  
  cat("Estimators:\n")
  cat("OLS\t", round(x$beta_ols, digits), "\n")
  if (!is.null(x$beta_tsls)) {
    cat("TSLS\t", round(x$beta_tsls, digits), "\n\n")
  } else {
    cat("\n")
  }
  
  cat(round((1-x$alpha)*100, digits), "% Confidence Intervals:\n", sep = "")
  cat("OLS\t[", round(x$confint_ols[1], digits), ",",
      round(x$confint_ols[2], digits), "]\n")
  if (!is.null(x$confint_tsls)) {
    cat("TSLS\t[", round(x$confint_tsls[1], digits), ",",
        round(x$confint_tsls[2], digits), "]\n\n")
  } else {
    cat("\n")
  }
  
  cat("Specified Bounds:\n")
  print(x$bounds, digits = digits)
  
  cat("\n")
  invisible(x)
}
