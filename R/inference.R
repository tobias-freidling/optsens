
one_opt <- function(data, indices, bounds, grid_specs, indep_x, dep_x,
                    print_warning, eps, extreme_value) {
  ## data as matrix with labelled columns
  ## create data_list
  y <- data[indices, "y"]
  d <- data[indices, "d"]
  z <- if ("z" %in% colnames(data)) data[indices, "z"] else NULL
  xp <- if (is.null(indep_x)) NULL else data[indices, indep_x]
  xt <- if (is.null(dep_x)) NULL else data[indices, dep_x]
  ##colnames(xp) <- indep_x
  ##colnames(xt) <- dep_x

  bounds_ind <- bootstrap_bounds(y, d, xt, xp, z, bounds, indices)

  grid_list <- feasible_grid(y, d, xt, xp, z, bounds_ind,
                             grid_specs, FALSE,
                             print_warning, eps)

  ## estimate ols and sds
  ylm_mat <- cbind(y, d, z, xp, xt)
  colnames(ylm_mat)[c(1,2)] <- c("y", "d")
  if (!is.null(z)) colnames(ylm_mat)[3] <- "z"

  dlm_mat <- cbind(d, z, xp, xt)
  colnames(dlm_mat)[1] <- "d"
  if (!is.null(z)) colnames(dlm_mat)[2] <- "z"

  ylm <- lm("y ~ -1 + .", data = as.data.frame(ylm_mat))
  dlm <- lm("d ~ -1 + .", data = as.data.frame(dlm_mat))
  beta_ols <- as.vector(ylm$coefficients["d"])
  sd_y_xzd <- sigma(ylm)
  sd_d_xz <- sigma(dlm)

  eval_mat <- eval_on_grid(grid_list$a_seq, grid_list$b_mat,
                           beta_ols, sd_y_xzd, sd_d_xz)

  all_na <- all(is.na(eval_mat))
  ret <- if(all_na) c(NA, NA) else c(min(eval_mat, na.rm = TRUE),
                                     max(eval_mat, na.rm = TRUE))
  if (!is.null(extreme_value)) {
    ret <- if (all_na) c(ret, (-1) * extreme_value, extreme_value) else c(ret, ret)
  }
  ret
}

#' @export
pir <- function(sa, grid_specs = list(num_x = 100, num_y = 100, num_z = 100),
                print_warning = FALSE, eps = 0.001) {
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
  pir <- one_opt(data_mat, indices, bounds, grid_specs, indep_x, dep_x,
                 print_warning, eps, NULL)
  pir
}

## other options: parallel = c("no", "multicore", "snow")
## ncpus

#' @export
sensint <- function(sa, alpha = 0.05, boot_procedure = c("perc", "basic"),
                    boot_samples = 500,
                    extreme_value = 1000,
                    grid_specs = list(num_x = 100, num_y = 100, num_z = 100),
                    print_warning = FALSE, eps = 0.001,
                    parallel = "no", ncpus = getOption("boot.ncpus", 1L)) {

  list2env(sa, environment())

  data_mat <- cbind(y, d, z, xp, xt)
  if (is.null(z)) {
    colnames(data_mat)[c(1,2)] <- c("y", "d")
  } else {
    colnames(data_mat)[c(1,2,3)] <- c("y", "d", "z")
  }
  indep_x <- if (is.null(xp)) NULL else colnames(xp)
  dep_x <- if (is.null(xt)) NULL else colnames(xt)

  boot_obj <- boot::boot(data = data_mat,
                         statistic = one_opt,
                         R = boot_samples,
                         bounds = bounds,
                         grid_specs = grid_specs,
                         indep_x = indep_x, dep_x = dep_x,
                         print_warning = print_warning,
                         eps = eps,
                         extreme_value = extreme_value,
                         parallel = parallel,
                         ncpus = ncpus)

  n_empty <- sum(is.na(boot_obj$t[,1]))

  lower <- boot::boot.ci(boot_obj, conf = 1 - alpha,
                         type = boot_procedure, index = 1)
  upper <- boot::boot.ci(boot_obj, conf = 1 - alpha,
                         type = boot_procedure, index = 2)
  lower_con <- boot::boot.ci(boot_obj, conf = 1 - alpha,
                             type = boot_procedure, index = 3)
  upper_con <- boot::boot.ci(boot_obj, conf = 1 - alpha,
                             type = boot_procedure, index = 4)

  ## Translation between boot_procedure and name in results
  boot_names <- c("normal", "basic", "student", "percent", "bca")
  names(boot_names) <- c("norm", "basic", "stud", "perc", "bca")

  sensint <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(sensint) <- c("sl", "su", "bootstrap", "conservative")

  ## extracting sensitivity intervals
  for (type in boot_names[boot_procedure]) {
    sl <- if (type == "normal") lower[[type]][2] else lower[[type]][4]
    su <- if (type == "normal") upper[[type]][3] else upper[[type]][5]
    sl_con <- if (type == "normal") lower_con[[type]][2] else lower_con[[type]][4]
    su_con <- if (type == "normal") upper_con[[type]][3] else upper_con[[type]][5]

    sensint[nrow(sensint) + 1, ] <- c(sl, su, type, FALSE)
    sensint[nrow(sensint) + 1, ] <- c(sl_con, su_con, type, TRUE)
  }


  list(boot_obj = boot_obj,
       n_empty = n_empty,
       sensint = sensint,
       alpha = alpha)
}
