
## Helper function creating the constraint functions
## Note: the abbreviations, e.g. a, b, d..., may not match the notation
## in the paper
constraint_functions <- function(y, d, xt, xp, z, bounds) {
  xtz <- cbind(xt, z)
  xtxp <- cbind(xt, xp)

  c3 <- r(y, d, cbind(xtz, xp))
  ## returns b
  eq_b <- function(a, g) (g - c3 * a) / sqrt(1-c3^2) / sqrt(1-a^2)
  ## returns g (c5, c6 are R^2)
  eq_g <- function(a, h, c4, c5, c6) {
    (h * sqrt(1-c4^2) * sqrt(1 - a^2 * (1-c5))
     + c4 * sqrt(1-c5) * a^2) / sqrt(1-c6)
  }

  ## Function to compute comparative-d bounds on g and comparative ZY-bounds
  c4 <- c5 <- c6 <- lb_h <- ub_h <- c()
  cyz <- bzy <- c()
  for (i in seq_len(nrow(bounds))) {
    if (bounds[i, "arrow"] == "UY" &&
        bounds[i, "kind"] == "comparative-d") {
      I <- unlist(bounds[i, "I"])
      xpic <- xp[, setdiff(colnames(xp), I)]
      xxi <- if (is.null(I)) xtz else cbind(xtz, xp[, I])
      c4 <- c(c4, r(y, d, xxi))
      c5 <- c(c5, r2(d, xpic, xxi))
      c6 <- c(c6, r2(y, xpic, xxi))
      lb_h <- c(lb_h, bounds[i, "lb"])
      ub_h <- c(ub_h, bounds[i, "ub"])
    } else if (bounds[i, "arrow"] == "ZY" &&
               bounds[i, "kind"] == "comparative") {
      J <- unlist(bounds[i, "J"])
      if (dim(xp)[2] < 2) {
        cyz <- c(cyz, r(y, xp, cbind(xtz, d))) ## only one xp -> J is disregarded
      } else {
        xpjc <- xp[, setdiff(colnames(xp), J)]
        cyz <- c(cyz, r(y, xp[, J], cbind(xtz, xpjc, d)))
      }
      bzy <- c(bzy, bounds[i, "b"])
    }
  }

  comp_d_bound_g <- function(a) {
    ub_g <- if (is.null(c4)) 1 else min(1, eq_g(a, ub_h, c4, c5, c6))
    lb_g <- if (is.null(c4)) -1 else max(-1, eq_g(a, lb_h, c4, c5, c6))
    c(lb_g, ub_g)
  }

  fun_list <- c(eq_b, comp_d_bound_g)
  names(fun_list) <- c("eq_b", "comp_d_bound_g")

  if (any(bounds[,"arrow"] %in% c("ZU", "ZY"))) {
    c1 <- r(y, z, cbind(xtxp, d))
    f_c1 <- c1 / sqrt(1-c1^2)
    c2 <- r(d, z, xtxp)
    ## returns fd
    eq_fd <- function(a, e) (e / sqrt(1-e^2) * sqrt(1-c2^2) - c2 * a) / sqrt(1-a^2)
    f_cyz <- if(is.null(cyz)) NULL else cyz / sqrt(1 - cyz^2)
    cyz_os <- if (is.null(cyz)) NULL else 1 - cyz^2

    ## compute bounds on f
    comp_bound_f <- function(a, b) {
      if (is.null(cyz)) {
        lb_f <- -1
        ub_f <- 1
      } else {
        f_comp <- (f_cyz * sqrt(1-a^2) + b * cyz * a) /
          sqrt(1-b^2) / sqrt(1 - a^2 * cyz_os)
        r2_comp <- f_comp^2 / (1 + f_comp^2)
        ub_f <- min(1, sqrt(bzy * r2_comp))
        lb_f <- max(-1, -sqrt(bzy * r2_comp))
      }
      c(lb_f, ub_f)
    }

    fun_list[[3]] <- eq_fd
    fun_list[[4]] <- f_c1
    fun_list[[5]] <- comp_bound_f
    names(fun_list)[3:5] <- c("eq_fd", "f_c1", "comp_bound_f")
  }
  fun_list
}



## Creating a grid a feasible (a,b)-values
## full_grid specifies if only values relevant for optimization, i.e.
## the boundary, are included.
#' @export
feasible_grid <- function(y, d, xt, xp, z, bounds, grid_specs, full_grid,
                          print_warning = FALSE, eps = 0.001) {
  list2env(grid_specs, environment())
  list2env(constraint_functions(y, d, xt, xp, z, bounds), environment())

  tsls_constr <- any(bounds[,"arrow"] %in% c("ZU", "ZY"))
  comp_uy_bound <- dim(bounds[(bounds$arrow == "UY") &
                                (bounds$kind != "direct"),])[1] > 0

  ## TSLS-only part
  if (tsls_constr) {
    ## Bound on e
    be <- bounds[bounds$arrow == "ZU", c("lb", "ub")]
    lb_e <- max(-1 + eps, be$lb)
    ub_e <- min(1 - eps, be$ub)

    ## Parameter-independent bounds on f (s for static)
    bf <- bounds[(bounds$arrow == "ZY") & (bounds$kind == "direct"), c("lb", "ub")]
    lb_fs <- max(-1 + eps, bf$lb)
    ub_fs <- min(1 - eps, bf$ub)
  }

  ## Bounds on a
  ba <- bounds[bounds$arrow == "UD", c("lb", "ub")]
  lb_a <- max(-1 + eps, ba$lb)
  ub_a <- min(1 - eps, ba$ub)

  ## Parameter-independent bounds on b and g (s for static)
  ## Bounds on b
  bb <- bounds[(bounds$arrow == "UY") & (bounds$kind == "direct"), c("lb", "ub")]
  lb_bs <- max(-1 + eps, bb$lb)
  ub_bs <- min(1 - eps, bb$ub)
  ## Bounds on g
  bg <- bounds[(bounds$arrow == "UY") & (bounds$kind == "comparative"), c("lb", "ub")]
  lb_gs <- max(-1 + eps, bg$lb)
  ub_gs <- min(1 - eps, bg$ub)

  a_seq <- seq(lb_a, ub_a, length.out = num_x)
  dim2 <- if (full_grid) num_y else 2
  b_mat <- matrix(NA, num_x, dim2)

  ## Grid-search
  for (i in seq_along(a_seq)) {
    ## Bound on b
    lb_b <- lb_bs
    ub_b <- ub_bs
    if (comp_uy_bound) {
      bg <- comp_d_bound_g(a_seq[i])
      lb_g <- max(lb_gs, bg[1])
      ub_g <- min(ub_gs, bg[2])
      lb_b <- max(lb_b, eq_b(a_seq[i], lb_g))
      ub_b <- min(ub_b, eq_b(a_seq[i], ub_g))
    }
    if (lb_b >= ub_b) {
      b_mat[i, ] <- NA
    } else if ((lb_b < ub_b) && (!tsls_constr)) {
      b_mat[i, ] <- if (full_grid) seq(lb_b, ub_b, length.out = num_y) else c(lb_b, ub_b)
    } else {
      ## case: tsls_constr
      fdl <- eq_fd(a_seq[i], lb_e)
      fdu <- eq_fd(a_seq[i], ub_e)
      lb_d <- fdl / sqrt(1 + fdl^2)
      ub_d <- fdu / sqrt(1 + fdu^2)
      d_seq <- seq(lb_d, ub_d, length.out = num_z)

      ## function that tests whether particular b is feasible
      test_b <- function (b) {
        bm <- comp_bound_f(a_seq[i], b)
        lb_f <- max(lb_fs, bm[1])
        ub_f <- min(ub_fs, bm[2])
        if (lb_f < ub_f) {
          ff <- (f_c1 * sqrt(1-d_seq^2) - b * d_seq) / sqrt(1-b^2)
          f <- ff / sqrt(1 + ff^2)
          ret <- any((f >= lb_f) & (f <= ub_f))
        } else {
          ret <- FALSE
        }
        ret
      }

      b_seq <- seq(lb_b, ub_b, length.out = num_y)
      if (full_grid) {
        for (j in seq_along(b_seq)) {
          b_mat[i, j] <- if(test_b(b_seq[j])) b_seq[j] else NA
        }
      } else {
        found <- FALSE
        for (j in seq_along(b_seq)) {
          found <- test_b(b_seq[j])
          if (found) {
            b_mat[i, 1] <- b_seq[j]
            break
          }
        }
        if (!found) {
          b_mat[i, ] <- NA
        } else {
          for (j in rev(seq_along(b_seq))) {
            if (test_b(b_seq[j])) {
              b_mat[i, 2] <- b_seq[j]
              break
            }
          }
        }
      }
    }
  }

  if (print_warning && all(is.na(b_mat))) {
    warning(paste0("Feasible values of R[D ~ U | X, Z] and ",
                   "R[Y ~ U | X, Z, D] could not be found.\n",
                   "The bounds are probably too restrictive."))
  }

  list(a_seq = a_seq, b_mat = b_mat)
}


## Evaluating objective, i.e. causal effect, over grid
eval_on_grid <- function(a_seq, b_mat, beta_ols, sd_y_xzd, sd_d_xz) {
  a_mat <- matrix(rep(a_seq, dim(b_mat)[2]), length(a_seq), dim(b_mat)[2])
  beta_ols - b_mat * a_mat / sqrt(1 - a_mat^2) * sd_y_xzd / sd_d_xz
}
