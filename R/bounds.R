
comp_bound_udy <- function(a, xp, xt, z, b, I, J) {
  xpic <- xp[, setdiff(colnames(xp), I)]
  xtpi <- if (is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])
  r2_axj_xxi <- if (dim(xp)[2] < 2) r2(a, xp, xtpi) else r2(a, xp[, J], xtpi)
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
        if (dim(xp)[2] < 2) {
          r2_yxj_xxI <- r2(y, xp, cbind(xtpi, d))
        } else {
          r2_yxj_xxI <- r2(y, xp[, J], cbind(xtpi, d))
        }
        bd <- sqrt(b * r2_yxj_xxI)
        b_max <- 1 / r2_yxj_xxI
        lb <- -bd
        ub <- bd
        if (bd >= 1 && print_warning) bound_warning(b_max)
      }
    } else if (arrow == "ZU") {
      if (dim(xp)[2] < 2) {
        r2_zxj_xxi <- r2(z, xp, xt)
      } else {
        xjc <- xp[, setdiff(colnames(xp), J)]
        r2_zxj_xxi <- r2(z, xp[, J], cbind(xt, xjc))
      }
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
  list(arrow = arrow,
       kind = kind,
       lb = lb, ub = ub, b = b,
       I = list(I),
       J = list(J))
}


## Adding bounds
## values of `arrow`: UD, UY, ZU, ZY
## values of `kind`: direct, comparative, comparative-d
#' @export
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


#' @export
remove_bound <- function(sa, name) {
  sa$bounds <- sa$bounds[rownames(sa$bounds) != name, ]
  sa
}


#' @export
view_bounds <- function(sa) {
  print(sa$bounds)
}


bootstrap_bounds <- function(y, d, xt, xp, z, bounds) {
  bounds_boot <- data.frame(matrix(nrow = 0, ncol = 7))
  colnames(bounds_boot) <- c("arrow", "kind", "lb", "ub", "b", "I", "J")
  for (i in seq_len(nrow(bounds))) {
    row <- compute_bound(y, d, xt, xp, z,
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
