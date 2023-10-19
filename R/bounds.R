## add_constraint
## remove_constraint
## view_constraints

## name, arrow, lb, ub, eta, b, I, J, partial_out_d
## arrow values: UD, UY, ZU, ZY
## kind values: direct, comparative, comparative-d

comp_bound_udy <- function(a, xp, xt, z, b, I, J) {
  xpic <- xp[, setdiff(colnames(xp), I)]
  xtpi <- if (is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])
  r2_axj_xxi <- r2(a, xp[, J], xtpi)
  r2_axic_xxi <- r2(a, xpic, xtpi)
  bd <- sqrt(b * r2_axj_xxi / (1 - r2_axic_xxi))
  b_max <- (1 - r2_axic_xxi) / r2_axj_xxi
  c(bd, b_max)
}

bound_warning <- function(b_max) {
  warning(paste0("The comparative bound must not exceed 1.\n",
                 "b must be smaller than ", b_max, "."))
}


compute_bound <- function(y, d, xt, xp, z,
                          arrow, kind, lb, ub, b, I, J, print_warning) {
  if (kind != "direct") {
    if (arrow == "UD") {
      res <- comp_bound_udy(d, xp, xt, z, b, I, J)
      lb <- -res[1]
      ub <- res[1]
      if (res[1] >= 1 && print_warning) bound_warning(res[2])
    } else if (arrow == "UY") {
      if (kind == "comparative") {
        res <- comp_bound_udy(y, xp, xt, z, b, I, J)
        lb <- -res[1]
        ub <- res[1]
        if (res[1] >= 1 && print_warning) bound_warning(res[2])
      } else {
        xtpi <- if(is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])
        r2_yxj_xxI <- r2(y, xp[, J], cbind(xtpi, d))
        bd <- sqrt(b * r2_yxj_xxI)
        b_max <- 1 / r2_yxj_xxI
        lb <- -bd
        ub <- bd
        if (bd >= 1 && print_warning) bound_warning(b_max)
      }
    } else if (arrow == "ZU") {
      xjc <- xp[, setdiff(colnames(xp), J)]
      r2_zxj_xxi <- r2(z, xp[, J], cbind(xt, xjc))
      bd <- b * r2_zxj_xxi
      b_max <- 1 / r2_zxj_xxi
      if(bd >= 1) {
        if(print_warning) bound_warning(b_max)
        lb <- -1 ## default bounds because actual expression may not be defined
        ub <- 1
      } else {
        b_zu_x <- sqrt(bd * (1 - r2_zxj_xxi) / (1 - r2_zxj_xxi * bd))
        lb <- -b_zu_x
        ub <- b_zu_x
      }
    }
  }

  ## J <- if (arrow == "ZY" && kind == "comparative") J else NA
  ## I <- if (arrow == "UY" && kind == "comparative-d") list(I) else NA
  list(arrow = arrow,
       kind = kind,
       lb = lb, ub = ub, b = b,
       I = list(I),
       J = list(J))
}


add_bound <- function(sa, arrow, kind,
                      lb = -1, ub = 1, b = NA, I = NULL, J = NULL,
                      name = "", print_warning = FALSE) {
  list2env(sa, environment())

  row_name <- if (name == "") paste0("b", rbc+1) else name
  row <- compute_bound(y, d, xt, xp, z,
                       arrow, kind, lb, ub, b, I, J, print_warning)
  sa$bounds[row_name, ] <- row
  sa$rbc <- rbc + 1
  sa
}


add_bound_old <- function(sa, arrow, kind,
                      lb = -1, ub = 1, b = NA, I = "", J = "",
                      correct_bound = FALSE, name = "", print_warning = FALSE) {
  list2env(sa, environment())
  if (arrow == "UD") {
    if (kind == "direct") {
      b <- NA
      J <- ""
    } else {
      xpic <- xp[, setdiff(colnames(xp), I)]
      xtpi <- if (is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])
      r2_dxj_xxi <- r2(d, xp[, J], xtpi)
      r2_dxic_xxi <- r2(d, xpic, xtpi)
      bd <- b * r2_dxj_xxi / (1 - r2_dxic_xxi)
      b_max <- (1 - r2_dxic_xxi) / r2_dxj_xxi

      if (bd >= 1) {
        if (correct_bound) {
          b <- b_max - 0.01
          bd <- b * r2_dxj_xxi / (1 - r2_dxic_xxi)
        } else if (print_warning) {
          warning(paste0("The comparative bound must not exceed 1.\n",
                         "For the given J and I, b must be ",
                         "smaller than ", b_max, "."))
        }
      }
      kind <- "comparative"
      lb <- -sqrt(bd)
      ub <- sqrt(bd)
      J <- ""
    }
  } else if (arrow == "UY") {
    if (kind == "direct") {
      b <- NA
      J <- ""
    } else if (kind == "comparative") {
      xpic <- xp[, setdiff(colnames(xp), I)]
      xtpi <- if(is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])

      r2_yxj_xxi <- r2(y, xp[, J], xtpi)
      r2_yxic_xxi <- r2(y, xpic, xtpi)
      bd <- b * r2_yxj_xxi / (1 - r2_yxic_xxi)
      b_max <- (1 - r2_yxic_xxi) / r2_yxj_xxi
      if (bd >= 1) {
        if (correct_bound) {
          b <- b_max - 0.01
          bd <- b * r2_yxj_xxi / (1 - r2_yxic_xxi)
        } else if (print_warning) {
          warning(paste0("The comparative bound must not exceed 1.\n",
                         "For the given J and I, b must be ",
                         "smaller than ", b_max, "."))
        }
      }
      lb <- -sqrt(bd)
      ub <- sqrt(bd)
      J <- ""
    } else {
      xtpi <- if(is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])
      r2_yxj_xxI <- r2(y, xp[, J], cbind(xtpi, d))
      bd <- b * r2_yxj_xxI
      b_max <- 1 / r2_yxj_xxI
      if (bd >= 1 && print_warning) {
        warning(paste0("The comparative bound must not exceed 1.\n",
                       "For the given J and I, b must be ",
                       "smaller than ", b_max, "."))
      }

      lb <- -sqrt(bd)
      ub <- sqrt(bd)
      J <- ""
      kind <- "comparative-d"
    }

  } else if (arrow == "ZU") {
    if(kind == "direct") {
      b <- NA
      J <- ""
    } else {
      xjc <- xp[, setdiff(indep_covariates, J)]
      r2_zxj_xxi <- r2(z, xp[, J], cbind(xt, xjc))
      bd <- b * r2_zxj_xxi
      b_max <- 1 / r2_zxj_xxi
      if(bd >= 1) {
        lb <- -1
        ub <- 1
        if(print_warning) {
          warning(paste0("The comparative bound must not exceed 1.\n",
                         "For the given J, b must be smaller ",
                         "than ", b_max, "."))
        }
      } else {
        b_zu_x <- sqrt(bd * (1 - r2_zxj_xxi) / (1 - r2_zxj_xxi * bd))
        lb <- -b_zu_x
        ub <- b_zu_x
      }
      J <- ""
      kind <- "comparative"
    }
  } else {
    if (kind == "direct") {
      b <- NA
      J <- ""
    } else {
      kind <- "comparative"
    }
    arrow <- "ZY"
  }


  row_name <- if (name == "") paste0("b", rbc) else name
  sa$bounds[row_name, ] <- list(arrow, kind, lb, ub, b, J)
  sa$rbc <- rbc + 1
  sa
}




remove_bound <- function(sa, name) {
  sa$bounds <- sa$bounds[rownames(sa$bounds) != name, ]
  sa
}


view_bounds <- function(sa) {
  print(sa$bounds)
}


bootstrap_bounds <- function(y, d, xt, xp, z, bounds, indices) {
  bounds_boot <- data.frame(matrix(nrow = 0, ncol = 7))
  colnames(bounds_boot) <- c("arrow", "kind", "lb", "ub", "b", "I", "J")

  for (i in seq_len(nrow(bounds))) {
    row <- compute_bound(y[indices], d[indices],
                         if (is.null(xt)) NULL else xt[indices,],
                         if (is.null(xp)) NULL else xp[indices,],
                         if (is.null(z)) NULL else z[indices],
                         bounds[i, "arrow"],
                         bounds[i, "kind"],
                         bounds[i, "lb"],
                         bounds[i, "ub"],
                         bounds[i, "b"],
                         unlist(bounds[i, "I"]),
                         unlist(bounds[i, "J"]),
                         FALSE)
    bounds_boot[paste0("b", i), ] <- row
  }

  bounds_boot
}
