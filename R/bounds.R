## Adding bounds
## values of `arrow`: UD, UY, ZU, ZY
## values of `kind`: direct, comparative, comparative-d
## depending on kind on arrow non-relevant parameters are ignored
#' @export
add_bound <- function(sa, arrow = c("UD", "UY", "ZU", "ZY"),
                      kind = c("direct", "comparative", "comparative-d"),
                      lb = -1, ub = 1, b = 1, I = NULL, J = NULL,
                      name = "", print_warning = TRUE) {
  
  arrow <- match.arg(arrow)
  kind <- match.arg(kind)
  
  check_add_bound(sa, arrow, kind, lb, ub, b, I, J, name, print_warning)
  ## Check all the inputs, e.g. "UD" + comparative-d doesn't work
  ## for IV bounds J must be a single item
  ## I and J must be disjoint
  ## I ignored for ZU, ZY
  
  ## warning that J will be ignored if dim(xp)[2] == 1
  
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
  
  if (class(sa) != "sensana") {
    stop("'sa' must be a 'sensana' object.")
  }
  
  if (!is.character(name) || length(name) != 1 ||
      !(name %in% rownames(sa$bounds))) {
    stop("'name' is not a string or not among the names of the specified bounds.")
  }
  
  sa$bounds <- sa$bounds[rownames(sa$bounds) != name, ]
  sa
}




comp_bound_udy <- function(a, xp, xt, z, b, I, J) {
  xpic <- xp[, setdiff(colnames(xp), I)]
  xtpi <- if (is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])
  xpj <- if(dim(xp)[2] < 2) xp else xp[,J]
  
  ##r2_axj_xxi <- if (dim(xp)[2] < 2) r2(a, xp, xtpi) else r2(a, xp[, J], xtpi)
  r2_axj_xxi <- r2(a, xpj, xtpi)
  r2_axic_xxi <- r2(a, xpic, xtpi)
  
  bd <- sqrt(b * r2_axj_xxi / (1 - r2_axic_xxi))
  b_max <- (1 - r2_axic_xxi) / r2_axj_xxi
  c(bd, b_max)
}



comp_d_bound_uy <- function(y, d, xp, xt, z, b, I, J) {
  xtpi <- if(is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])
  xpj <- if(dim(xp)[2] < 2) xp else xp[, J]
  r2_yxj_xxI <- r2(y, xpj, cbind(xtpi, d))
  
  bd <- sqrt(b * r2_yxj_xxI)
  b_max <- 1 / r2_yxj_xxI
  c(bd, b_max)
}



comp_bound_zu <- function(xp, xt, z, b, J) {
  
  if (dim(xp)[2] < 2) {
    r2_zxj_xxi <- r2(z, xp, xt)
  } else {
    r2_zxj_xxi <- r2(z, xp[, J], cbind(xt, xp[, setdiff(colnames(xp), J)]))
  }
  
  bd <- b * r2_zxj_xxi
  b_max <- 1 / r2_zxj_xxi
  
  if (b >= b_max) {
    ## default bounds because actual expression may not be defined
    res <- c(1, b_max)
  } else {
    res <- c(sqrt(bd * (1 - r2_zxj_xxi) / (1 - r2_zxj_xxi * bd)), b_max)
  }
  
  res
}



bound_warning <- function(b_max) {
  warning(paste0("The comparative bound exceeds 1.\n",
                 "For a bound smaller than 1, b must be smaller than ", b_max, "."))
}



compute_bound <- function(y, d, xt, xp, z,
                           arrow, kind, lb, ub, b, I, J, print_warning) {
  
  if (kind == "direct") {
    bounds <- c(lb, ub)
    b <- NA
    I <- J <- NULL
  } else {
    if (arrow == "UD") {
      res <- comp_bound_udy(d, xp, xt, z, b, I, J)
    } else if (arrow == "UY" && kind == "comparative") {
      res <- comp_bound_udy(y, xp, xt, z, b, I, J)
    } else if (arrow == "UY" && kind == "comparative-d") {
      res <- comp_d_bound_uy(y, d, xp, xt, z, b, I, J)
    } else if (arrow == "ZU") {
      res <- comp_bound_zu(xp, xt, z, b, J)
      I <- NULL
    } else { ## arrow == "ZY"
      res <- NA
      I <- NULL
    }
    bounds <- c(-res[1], res[1])
    if (!is.na(res[1]) && res[1] >= 1 && print_warning) bound_warning(res[2])
  }
  
  list(arrow = arrow, kind = kind, lb = bounds[1], ub = bounds[2],
       b = b, I = list(I), J = list(J))
}



recompute_bounds <- function(y, d, xt, xp, z, bounds) {
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
