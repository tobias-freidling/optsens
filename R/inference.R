
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


## Is there a smarter way to handle the data hand-over??

## Solving the optimization problem with the observations data[indices,]
one_opt <- function(data, indices, bounds, grid_specs, indep_x, dep_x, eps, cpp) {
  ## data as matrix with labelled columns
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

  ## Core of the method
  bounds_ind <- recompute_bounds(y, d, xt, xp, z, bounds)
  grid_list <- feasible_grid(y, d, xt, xp, z, bounds_ind,
                             grid_specs, FALSE, eps, cpp)
  eval_param <- compute_eval_param(y, d, xt, xp, z)
  eval_mat <- eval_on_grid(grid_list$a_seq, grid_list$b_mat, eval_param)

  ## Evaluating and reporting the results
  ret <- if(all(is.na(eval_mat))) c(-Inf, Inf) else range(eval_mat, na.rm = TRUE)
  ret
}


## Computing the partially identified region (PIR)
#' @export
pir <- function(sa, grid_specs = list(num_x = 100, num_y = 100, num_z = 100),
                eps = 0.001, cpp = FALSE) {
  
  check_pir(sa, grid_specs, eps)
  list2env(sa, environment())

  data_mat <- cbind(y, d, z, xp, xt)
  if (is.null(z)) {
    colnames(data_mat)[c(1,2)] <- c("y", "d")
  } else {
    colnames(data_mat)[c(1,2,3)] <- c("y", "d", "z")
  }

  indices <- 1:length(y)
  indep_x <- if (is.null(xp)) NULL else colnames(xp)
  dep_x <- if (is.null(xt)) NULL else colnames(xt)
  pir <- one_opt(data_mat, indices, bounds, grid_specs, indep_x, dep_x, eps, cpp)
  pir
}





## Computing a (1-alpha) sensitivity interval via the bootstrap
#' @export
sensint <- function(sa, alpha = 0.05, boot_procedure = c("perc", "basic"),
                    boot_samples = 500,
                    grid_specs = list(num_x = 100, num_y = 100, num_z = 100),
                    eps = 0.001,
                    parallel = "no", ncpus = getOption("boot.ncpus", 1L)) {

  check_sensint(sa, alpha, boot_procedure, boot_samples,
                grid_specs, eps, parallel, ncpus)
  
  list2env(sa, environment())
  
  
  ## boot_procedure
  
  ## warning for bca
  
  ## how to do match.arg with possibly vector-valued input?

  data_mat <- cbind(y, d, z, xp, xt)
  if (is.null(z)) {
    colnames(data_mat)[c(1,2)] <- c("y", "d")
  } else {
    colnames(data_mat)[c(1,2,3)] <- c("y", "d", "z")
  }
  indep_x <- if (is.null(xp)) NULL else colnames(xp)
  dep_x <- if (is.null(xt)) NULL else colnames(xt)

  ## Bootstrapping the one_opt function
  boot_obj <- boot::boot(data = data_mat,
                         statistic = one_opt,
                         R = boot_samples,
                         bounds = bounds,
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



#' @export
print.sensint <- function(sint_obj, digits = max(3L, getOption("digits") - 3L)) {
  list2env(sint_obj, environment())
  
  cat(round((1-alpha)*100, digits), "% Sensitivity Intervals:\n", sep = "")
  print(sensint, digits)
  cat("\n")
  
  invisible(sint_obj)
}
