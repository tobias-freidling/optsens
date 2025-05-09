
## Helper function computing parameters to evaluate the objective
compute_eval_param <- function(y, d, xt, xp, z) {
  cov_mat <- cbind(xt, xp, z)
  ylm <- stats::lm.fit(cbind(d, cov_mat), y)
  dlm <- stats::lm.fit(cov_mat, d)

  beta_ols <- as.vector(ylm$coefficients[1])
  sd_y_xzd <- sqrt(sum(ylm$residuals^2) / (dim(cov_mat)[1] - ylm$rank))
  sd_d_xz <- sqrt(sum(dlm$residuals^2) / (dim(cov_mat)[1] - dlm$rank))

  list(beta_ols = beta_ols, sd_y_xzd = sd_y_xzd, sd_d_xz = sd_d_xz)
}



## Solving the optimization problem with the observations data[indices,]
## written in this way to facilitate using the boot-package
one_opt <- function(data, indices, bounds, grid_specs, indep_x, dep_x, eps) {
  ## Disentangling data to pass it on to the core methods
  ## (data is a matrix with labelled columns)
  y <- data[indices, "y"]
  d <- data[indices, "d"]
  z <- if ("z" %in% colnames(data)) data[indices, "z"] else NULL

  if (is.null(indep_x)) {
    xp <- NULL
  } else if (length(indep_x) == 1) {
    xp <- as.matrix(data[indices, indep_x])
    colnames(xp) <- indep_x
  } else {
    xp <- data[indices, indep_x]
  }

  if (is.null(dep_x)) {
    xt <- NULL
  } else if (length(dep_x) == 1) {
    xt <- as.matrix(data[indices, dep_x])
    colnames(xt) <- dep_x
  } else {
    xt <- data[indices, dep_x]
  }

  ## Recomputing the bounds for the new indices
  bounds_ind <- recompute_bounds(y, d, xt, xp, z, bounds)
  ## Computing a grid of feasible points
  grid_list <- feasible_grid(y, d, xt, xp, z, bounds_ind,
                             grid_specs, FALSE, eps)
  
  ## Evaluation of the objective/causal effect over the grid
  ep <- compute_eval_param(y, d, xt, xp, z)
  eval_mat <- eval_on_grid(grid_list$p1_seq, grid_list$p2_mat,
                           ep$beta_ols, ep$sd_y_xzd, ep$sd_d_xz)
  
  ## Computing the smallest and largest value
  ret <- if(all(is.na(eval_mat))) c(-Inf, Inf) else range(eval_mat, na.rm = TRUE)
  ret
}


#' Compute the Partially Identified Range (PIR)
#'
#' Computes the PIR by finding the smallest and largest possible value of the
#' causal effect according to the specified sensitivity model. To solve this
#' optimization problem a tailored grid search algorithm is used.
#' 
#'
#' @param sa An object of the class \code{sensana}.
#' @param grid_specs A named list of the three numeric values `N1`, `N2` and `N5` specifying
#'   the number of points considered in the grid search for each of the three dimensions.
#' @param eps A small numeric number to bound sensitivity parameters away from `1`
#'   and `-1`. Default is \code{0.001}.
#'
#' @return A numeric vector of length `2`, the estimated partially identified range.
#'
#' @examples
#' set.seed(123)
#' u <- rnorm(20)
#' xt <- -0.25*u + rnorm(20)
#' xp <- 0.5 * xt + rnorm(20)
#' z <- 0.25 * xt - xp + u + rnorm(20)
#' d <- 0.5 * xt + 0.5 * xp + u + 2*z + rnorm(20)
#' y <- d + 2*xp - xt + z + u + rnorm(20)
#' sa <- sensana(y = y, d = d, x = data.frame(xt = xt, xp = xp), z = z,
#'               dep_x = "xt", indep_x = "xp", alpha = 0.05, quantile = "t")
#'
#' sa <- add_bound(sa, arrow = "UD", kind = "direct", lb = -0.3, ub = 0.3)
#' sa <- add_bound(sa, arrow = "ZY", kind = "direct", lb = -0.3, ub = 0.3)
#'                 
#' pir_res <- pir(sa, grid_specs = list(N1 = 100, N2 = 100, N5 = 100))
#'
#' @export
pir <- function(sa, grid_specs = list(N1 = 100, N2 = 100, N5 = 100),
                eps = 0.001) {
  
  check_pir(sa, grid_specs, eps)
  
  data_mat <- cbind(sa$y, sa$d, sa$z, sa$xp, sa$xt)
  if (is.null(sa$z)) {
    colnames(data_mat)[c(1,2)] <- c("y", "d")
  } else {
    colnames(data_mat)[c(1,2,3)] <- c("y", "d", "z")
  }

  indices <- 1:length(sa$y)
  indep_x <- if (is.null(sa$xp)) NULL else colnames(sa$xp)
  dep_x <- if (is.null(sa$xt)) NULL else colnames(sa$xt)
  pir <- one_opt(data_mat, indices, sa$bounds, grid_specs, indep_x, dep_x, eps)
  pir
}







#' Compute Sensitivity Intervals
#' 
#' @description 
#' Computes a `1-alpha` sensitivity interval for the partially identified range
#' via the bootstrap. Users can specify the bootstrap method(s) used for constructing
#' the sensitivity interval.
#' 
#' 
#' @details 
#' Since the sensitivity model given by the user-specified bounds is data-dependent,
#' it needs to be re-computed on every bootstrap sample. In cases where the resulting
#' sensitivity model is empty, we can either discard such a bootstrap sample or set
#' the estimated PIR to \eqn{(-\infty, \infty)}. In the returned data frame containing
#' the sensitivity intervals, these 2 possibilities are indicated via `conservative = FALSE`
#' and `conservative = TRUE`. We recommend using conservative sensitivity intervals.
#'
#' @param sa An object of the class \code{sensana}.
#' @param alpha Significance level used to construct the sensitivity interval.
#'   Default is \code{0.05} for a \eqn{95\%} sensitivity interval.
#' @param boot_procedure A vector of or a single character string specifying
#'   the bootstrap method(s) to construct the sensitivity interval. Accepted
#'   values are \code{"perc"} (percentile), \code{"basic"} (basic), \code{"bca"}
#'   (BCa), \code{"norm"} (normal approximation) and \code{"stud"} (studentized).
#' @param boot_samples Integer; the number of bootstrap resamples. Default is \code{1000}.
#' @param grid_specs A named list of the three numeric values `N1`, `N2` and `N5` specifying
#'   the number of points considered in the grid search for each of the three dimensions.
#' @param eps A small numeric number to bound sensitivity parameters away from `1`
#'   and `-1`. Default is \code{0.001}.
#' @param parallel A character string passed to \code{\link[boot]{boot}}
#'   specifying the parallel processing backend. Options are \code{"no"},
#'   \code{"multicore"}, or \code{"snow"}. Default is `"no"`.
#' @param ncpus Integer number of CPU cores to use for parallel processing.
#'   Passed to \code{\link[boot]{boot}}.   
#'   
#'   
#'   
#' @returns Object of the class `sensint`. This is a list
#'  containing the following components:
#'  \describe{
#'  \item{sensint}{A data frame where each row contains one sensitivity interval.
#'    The columns `sl` and `su` refer to the lower and upper end of the respective
#'    sensitivity interval. The column `bootstrap` specifies the method used to
#'    construct the interval. The column `conservative` indicates whether the
#'    sensitivity interval is conservative or not (see details above).}
#'  \item{alpha}{The significance level.}
#'  \item{n_empty}{The number of bootstrap resamples where the sensitivity model
#'    is empty.}
#'  \item{boot_obj}{The `boot` object returned by the internal call to [boot::boot()].}
#'  }
#' 
#'
#' @examples
#' set.seed(123)
#' u <- rnorm(20)
#' xt <- -0.25*u + rnorm(20)
#' xp <- 0.5 * xt + rnorm(20)
#' z <- 0.25 * xt - xp + u + rnorm(20)
#' d <- 0.5 * xt + 0.5 * xp + u + 2*z + rnorm(20)
#' y <- d + 2*xp - xt + z + u + rnorm(20)
#' sa <- sensana(y = y, d = d, x = data.frame(xt = xt, xp = xp), z = z,
#'               dep_x = "xt", indep_x = "xp", alpha = 0.05, quantile = "t")
#'
#' sa <- add_bound(sa, arrow = "UD", kind = "direct", lb = -0.3, ub = 0.3)
#' sa <- add_bound(sa, arrow = "ZY", kind = "direct", lb = -0.3, ub = 0.3)
#'                 
#' si <- sensint(sa, alpha = 0.05, boot_procedure = c("perc", "bca"),
#'               boot_samples = 5000)
#' 
#' @export
sensint <- function(sa, alpha = 0.05, boot_procedure = c("perc", "basic", "bca"),
                    boot_samples = 1000,
                    grid_specs = list(N1 = 100, N2 = 100, N5 = 100),
                    eps = 0.001,
                    parallel = "no", ncpus = getOption("boot.ncpus", 1L)) {

  check_sensint(sa, alpha, boot_procedure, boot_samples,
                grid_specs, eps, parallel, ncpus)

  data_mat <- cbind(sa$y, sa$d, sa$z, sa$xp, sa$xt)
  if (is.null(sa$z)) {
    colnames(data_mat)[c(1,2)] <- c("y", "d")
  } else {
    colnames(data_mat)[c(1,2,3)] <- c("y", "d", "z")
  }
  indep_x <- if (is.null(sa$xp)) NULL else colnames(sa$xp)
  dep_x <- if (is.null(sa$xt)) NULL else colnames(sa$xt)

  ## Bootstrapping the one_opt function
  boot_obj <- boot::boot(data = data_mat,
                         statistic = one_opt,
                         R = boot_samples,
                         bounds = sa$bounds,
                         grid_specs = grid_specs,
                         indep_x = indep_x,
                         dep_x = dep_x,
                         eps = eps,
                         parallel = parallel,
                         ncpus = ncpus)

  ## Reporting results

  n_empty <- sum(boot_obj$t[,2] == Inf)
  lam <- n_empty / boot_samples
  
  boot_names <- c("normal", "basic", "student", "percent", "bca")
  names(boot_names) <- c("norm", "basic", "stud", "perc", "bca")
  
  if (lam < alpha) {
    lower <- boot::boot.ci(boot_obj,
                           conf = c(1-alpha, 1-alpha+lam),
                           type = boot_procedure,
                           index = 1)
    
    upper <- boot::boot.ci(boot_obj,
                           conf = c(1-alpha, 1-alpha+lam),
                           type = boot_procedure,
                           index = 2)
    
    ## extracting sensitivity intervals
    sl <- su <- bootstrap <- conservative <- c()
    for (type in boot_names[boot_procedure]) {
      sl <- c(sl, if(type == "normal") lower[[type]][1:2,2] else lower[[type]][1:2,4])
      su <- c(su, if(type == "normal") upper[[type]][1:2,3] else upper[[type]][1:2,5])
      bootstrap <- c(bootstrap, rep(type,2))
      conservative <- c(conservative, c(FALSE, TRUE))
    }
  } else {
    lower <- boot::boot.ci(boot_obj,
                           conf = 1-alpha,
                           type = boot_procedure,
                           index = 1)
    
    upper <- boot::boot.ci(boot_obj,
                           conf = 1-alpha,
                           type = boot_procedure,
                           index = 2)
    
    sl <- su <- bootstrap <- conservative <- c()
    for (type in boot_names[boot_procedure]) {
      sl <- c(sl, if(type == "normal") lower[[type]][2] else lower[[type]][4])
      su <- c(su, if(type == "normal") upper[[type]][3] else upper[[type]][5])
      bootstrap <- c(bootstrap, type)
      conservative <- c(conservative, FALSE)
    }
  }
  
  sensint <- data.frame(sl = sl, su = su,
                        bootstrap = bootstrap,
                        conservative = conservative)
  ret <- list(sensint = sensint, alpha = alpha,
              n_empty = n_empty, boot_obj = boot_obj)
  class(ret) <- "sensint"
  ret
}



#' Print Method for `sensint` objects
#'
#' @param x An object of class \code{"sensint"}.
#' @param digits The number of significant digits to display.
#' @param ... Additional arguments passed to or from other methods (ignored here).
#'
#' @method print sensint
#' @return Invisibly returns \code{x}, the original object.
#' 
#' @examples
#' set.seed(123)
#' u <- rnorm(20)
#' xt <- -0.25*u + rnorm(20)
#' xp <- 0.5 * xt + rnorm(20)
#' z <- 0.25 * xt - xp + u + rnorm(20)
#' d <- 0.5 * xt + 0.5 * xp + u + 2*z + rnorm(20)
#' y <- d + 2*xp - xt + z + u + rnorm(20)
#' sa <- sensana(y = y, d = d, x = data.frame(xt = xt, xp = xp), z = z,
#'               dep_x = "xt", indep_x = "xp", alpha = 0.05, quantile = "t")
#'
#' sa <- add_bound(sa, arrow = "UD", kind = "direct", lb = -0.3, ub = 0.3)
#' sa <- add_bound(sa, arrow = "ZY", kind = "direct", lb = -0.3, ub = 0.3)
#'                 
#' sensint_obj <- sensint(sa, alpha = 0.05)
#' print(sensint_obj)
#' 
#' @export
print.sensint <- function(x, digits = max(3L, getOption("digits") - 3L),...) {
  cat(round((1-x$alpha)*100, digits), "% Sensitivity Intervals:\n", sep = "")
  print(x$sensint, digits)
  cat("\n")
  
  invisible(x)
}
