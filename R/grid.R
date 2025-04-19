

ols_constraint_functions <- function(y, d, xt, xp, z, bounds) {
  ## hb
  xtz <- cbind(xt, z)
  c1 <- r(y, d, cbind(xtz, xp))
  hb <- function(a, dpar) (dpar - c1 * a) / sqrt(1-c1^2) / sqrt(1-a^2)
  ## dpar instead of d to avoid ambiguity
  
  ## hd
  hd <- function(a, e, c2, c3, c4) {
    (e * sqrt(1-c2^2) * sqrt(1 - a^2 * (1-c3^2))
     + c2 * sqrt(1-c3^2) * a^2) / sqrt(1-c4^2)
  }
  
  
  ## Computing constants
  c2 <- c3 <- c4 <- lb_e <- ub_e <- c()
  for (i in seq_len(nrow(bounds))) {
    if (bounds[i, "arrow"] == "UY" &&
        bounds[i, "kind"] == "comparative-d") {
      I <- unlist(bounds[i, "I"])
      xpic <- xp[, setdiff(colnames(xp), I)]
      xxi <- if (is.null(I)) xtz else cbind(xtz, xp[, I])
      c2 <- c(c2, r(y, d, xxi))
      c3 <- c(c3, r(d, xpic, xxi))
      c4 <- c(c4, r(y, xpic, xxi))
      lb_e <- c(lb_e, bounds[i, "lb"])
      ub_e <- c(ub_e, bounds[i, "ub"])
    }
  }
  
  ## Function yielding comparative-d bounds on d
  comp_d_bound_d <- function(a) {
    if (is.null(c2)) {
      ret <- c(-1,1)
    } else {
      ret <- c(max(-1, hd(a, lb_e, c2, c3, c4)),
               min(1, hd(a, ub_e, c2, c3, c4)))
    }
    ret
  }
  
  fun_list <- c(hb, comp_d_bound_d)
  names(fun_list) <- c("hb", "comp_d_bound_d")
  fun_list
}


tsls_constraint_functions <- function(y, d, xt, xp, z, bounds) {
  xtz <- cbind(xt, z)
  x <- cbind(xt, xp)
  
  ## hfg
  c5 <- r(d, z, x)
  hfg <- function(a, m) (f(m) * sqrt(1-c5^2) - c5 * a) / sqrt(1-a^2)
  
  
  ## hfo
  c6 <- r(y, z, cbind(x, d))
  hfo <- function(b, g) (f(c6) * sqrt(1-g^2) - b * g) / sqrt(1-b^2)
  
  
  ## hfq
  hfq <- function(a, b, c7) {
    (f(c7) * sqrt(1-a^2) + c7 * a * b) / sqrt(1-b^2) / sqrt(1 - a^2 * (1-c7^2))
  }
  
  
  ## Computing constants
  c7 <- b7 <- c()
  for (i in seq_len(nrow(bounds))) {
    if (bounds[i, "arrow"] == "ZY" &&
        bounds[i, "kind"] == "comparative") {
      J <- unlist(bounds[i, "J"])
      if (dim(xp)[2] < 2) {
        c7 <- c(c7, r(y, xp, cbind(xtz, d))) ## only one xp -> J is disregarded
      } else {
        xpjc <- xp[, setdiff(colnames(xp), J)]
        c7 <- c(c7, r(y, xp[, J], cbind(xtz, xpjc, d)))
      }
      b7 <- c(b7, bounds[i, "b"])
    }
  }
  
  ## Function yielding comparative bounds on o
  comp_bound_o <- function(a, b) {
    if (is.null(c7)) {
      lb_o <- -1
      ub_o <- 1
    } else {
      fq <- hfq(a, b, c7)
      bounds_o <- sqrt(b7 * finv(fq)^2)
      ub_o <- min(1, bounds_o)
      lb_o <- max(-1, -bounds_o)
    }
    c(lb_o, ub_o)
  }
  
  fun_list <- list(hfg, hfo, comp_bound_o)
  names(fun_list) <- c("hfg", "hfo", "comp_bound_o")
  fun_list
}


constraint_functions <- function(y, d, xt, xp, z, bounds) {
  ## Functions for OLS constraints
  list_ols <- ols_constraint_functions(y, d, xt, xp, z, bounds)
  
  ## Functions for TSLS constraints
  if (any(bounds[,"arrow"] %in% c("ZU", "ZY"))) {
    list_tsls <- tsls_constraint_functions(y, d, xt, xp, z, bounds)
  } else {
    list_tsls <- NULL
  }
  
  c(list_ols, list_tsls)
}


static_bounds <- function(y, d, xt, xp, z, bounds, eps) {
  
  ## Bounds on a
  ba <- bounds[bounds$arrow == "UD", c("lb", "ub")]
  lb_as <- max(-1 + eps, ba$lb)
  ub_as <- min(1 - eps, ba$ub)
  
  ## Bounds on b
  bb <- bounds[(bounds$arrow == "UY") & (bounds$kind == "direct"), c("lb", "ub")]
  lb_bs <- max(-1 + eps, bb$lb)
  ub_bs <- min(1 - eps, bb$ub)
  
  ## Bounds on d
  bd <- bounds[(bounds$arrow == "UY") & (bounds$kind == "comparative"), c("lb", "ub")]
  lb_ds <- max(-1 + eps, bd$lb)
  ub_ds <- min(1 - eps, bd$ub)
  
  ## for TSLS constraints
  if (any(bounds[,"arrow"] %in% c("ZU", "ZY"))) {
    ## Bound on m
    bm <- bounds[bounds$arrow == "ZU", c("lb", "ub")]
    lb_ms <- max(-1 + eps, bm$lb)
    ub_ms <- min(1 - eps, bm$ub)
    
    ## Bounds on o
    bo <- bounds[(bounds$arrow == "ZY") & (bounds$kind == "direct"), c("lb", "ub")]
    lb_os <- max(-1 + eps, bo$lb)
    ub_os <- min(1 - eps, bo$ub)
  } else {
    lb_ms <- ub_ms <- lb_os <- ub_os <- NULL
  }
  
  list(lb_ms = lb_ms, ub_ms = ub_ms, lb_os = lb_os, ub_os = ub_os,
       lb_as = lb_as, ub_as = ub_as, lb_bs = lb_bs, ub_bs = ub_bs,
       lb_ds = lb_ds, ub_ds = ub_ds)
}


grid_search <- function(constraint_functions, static_bounds, grid_specs,
                        full_grid, exist_tsls_bound, exist_comp_uy_bound) {
  
  ## needs comparative bounds and TSLS info
  ## STOPPED HERE (good effort)
  
  list2env(constraint_functions, environment())
  list2env(static_bounds, environment())
  list2env(grid_specs, environment())
  
  
  a_seq <- seq(lb_as, ub_as, length.out = num_x)
  b_mat <- matrix(NA, num_x, if (full_grid) num_y else 2)
  
  ## Grid-search
  for (i in seq_along(a_seq)) {
    ## Push forward OLS bounds on b
    lb_b <- lb_bs
    ub_b <- ub_bs
    if (exist_comp_uy_bound) {
      bd <- comp_d_bound_d(a_seq[i])
      lb_d <- max(lb_ds, bd[1])
      ub_d <- min(ub_ds, bd[2])
      lb_b <- max(lb_b, hb(a_seq[i], lb_d))
      ub_b <- min(ub_b, hb(a_seq[i], ub_d))
    }
    
    if (lb_b >= ub_b) {
      b_mat[i, ] <- NA
    } else if ((lb_b < ub_b) && (!exist_tsls_bound)) {
      b_mat[i, ] <- if (full_grid) seq(lb_b, ub_b, length.out = num_y) else c(lb_b, ub_b)
    } else {
      ## Case with TSLS constraints
      
      ## Push forward ZU bounds on g
      lb_g <- finv(hfg(a_seq[i], lb_ms))
      ub_g <- finv(hfg(a_seq[i], ub_ms))
      g_seq <- seq(lb_g, ub_g, length.out = num_z)
      
      ## function that tests whether particular b is feasible
      test_b <- function (b) {
        bo <- comp_bound_o(a_seq[i], b)
        lb_o <- max(lb_os, bo[1])
        ub_o <- min(ub_os, bo[2])
        if (lb_o < ub_o) {
          o <- finv(hfo(b, g_seq))
          ret <- any((o >= lb_o) & (o <= ub_o))
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
  
  
  # if (print_warning && all(is.na(b_mat))) {
  #   warning(paste0("Feasible values of R[D ~ U | X, Z] and ",
  #                  "R[Y ~ U | X, Z, D] could not be found.\n",
  #                  "The bounds are probably too restrictive."))
  # }
  
  list(a_seq = a_seq, b_mat = b_mat)
}





feasible_grid <- function(y, d, xt, xp, z, bounds, grid_specs,
                          full_grid, eps = 0.001) {
  
  constraint_functions <- constraint_functions(y, d, xt, xp, z, bounds)
  static_bounds <- static_bounds(y, d, xt, xp, z, bounds, eps)
  
  exist_tsls_bound <- any(bounds[,"arrow"] %in% c("ZU", "ZY"))
  exist_comp_uy_bound <- dim(bounds[(bounds$arrow == "UY") &
                                      (bounds$kind != "direct"),])[1] > 0
  
  res <- grid_search(constraint_functions, static_bounds, grid_specs,
                     full_grid, exist_tsls_bound, exist_comp_uy_bound)
  res
}





## Evaluating objective, i.e. causal effect, over grid
eval_on_grid <- function(a_seq, b_mat, eval_param) {
  list2env(eval_param, environment())
  a_mat <- matrix(rep(a_seq, dim(b_mat)[2]), length(a_seq), dim(b_mat)[2])
  beta_ols - b_mat * f(a_mat) * sd_y_xzd / sd_d_xz
}
