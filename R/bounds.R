#' Add a Sensitivity Bound to a `"sensana"` object
#'
#' This method adds a new constraint / bound to the sensitivity model specified
#' in `sa`. Users can specify bounds on the effects of \eqn{U} on \eqn{D} (`"UD"`),
#' \eqn{U} on \eqn{Y} (`"UY"`), \eqn{Z} on \eqn{U} (`"ZU"`) or
#' \eqn{Z} on \eqn{Y} (`"ZY"`). There are 3 different availble kinds of bounds:
#' `"direct"`, `"comparative"` and `"comparative-d"`.
#' (The latter is only suitable for `"UY"`.)
#' 
#' 
#' @details
#' Arguments not pertaining to the specified `arrow` and `kind` are ignored.
#' If a comparative bound on `"ZU"` or `"ZY"` is given, `J` must be a single
#' string; moreover, it will be ignored if `sa$indep_x` contains only one entry.
#' 
#' For the interpretation of the bounds, see the article cited below.
#' 
#'
#' @param sa Object of the class \code{sensana}.
#' @param arrow String indicating the type of arrow to add.
#'   One of \code{"UD"}, \code{"UY"}, \code{"ZU"}, or \code{"ZY"}.
#' @param kind String specifying the type of the bound.
#'   One of \code{"direct"}, \code{"comparative"}, or \code{"comparative-d"}.
#' @param lb Numeric value specifying the lower bound. Default is \code{-1}.
#' @param ub Numeric value specifying the upper bound. Default is \code{1}.
#' @param b Positive numeric value used for comparative bounds. Default is \code{1}.
#' @param I Character vector of covariate indices/names used for comparative bounds;
#'   subset of `sa$indep_x`. Default is \code{NULL}.
#' @param J Character vector of covariate indices/names used for comparative bounds;
#'   subset of `sa$indep_x`. Default is \code{NULL}.
#' @param name Character string giving a name to the bound. Default is an empty string.
#' @param print_warning Logical; if \code{TRUE}, a warning will be printed if the bound
#'  is does not have an effect on the estimation of the partially identified range.
#'  Default is \code{TRUE}.
#'
#' @return An updated \code{sensana} object that includes the added bound in its
#'  `bounds` data frame.
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
#' ## Direct bounds
#' sa <- add_bound(sa, arrow = "UD", kind = "direct", lb = -0.3, ub = 0.3)
#' sa <- add_bound(sa, arrow = "ZY", kind = "direct", lb = -0.3, ub = 0.3)
#' 
#' ## Comparative bounds
#' sa <- add_bound(sa, arrow = "UD", kind = "comparative",
#'                 b = 0.2, J = "xp", I = NULL)
#' sa <- add_bound(sa, arrow = "UY", kind = "comparative-d",
#'                 b = 0.3, J = "xp", I = NULL)
#' sa <- add_bound(sa, arrow = "ZU", kind = "comparative",
#'                 b = 0.4, J = "xp")
#'
#' @references
#' Freidling T, Zhao Q (2024). “Optimization-based Sensitivity Analysis for
#' Unmeasured Confounding using Partial
#' Correlations.” _arXiv preprint arXiv:2301.00040v3_.
#' @export
add_bound <- function(sa, arrow = c("UD", "UY", "ZU", "ZY"),
                      kind = c("direct", "comparative", "comparative-d"),
                      lb = -1, ub = 1, b = 1, I = NULL, J = NULL,
                      name = "", print_warning = TRUE) {
  
  arrow <- match.arg(arrow)
  kind <- match.arg(kind)
  
  check_add_bound(sa, arrow, kind, lb, ub, b, I, J, name, print_warning)
  row_name <- if (name == "") paste0("b", sa$rbc+1) else name
  row <- compute_bound(sa$y, sa$d, sa$xt, sa$xp, sa$z,
                       arrow, kind, lb, ub, b, I, J, print_warning)
  sa$bounds[row_name, ] <- row
  sa$rbc <- sa$rbc + 1
  sa
}



#' Remove a Sensitivity Bound from a `"sensana"` object
#'
#' Removes a previously added sensitivity bound from a \code{sensana} object by name.
#' 
#' @details
#' To display the names of all the current bounds, print the `sensana` object or
#' directly access the `sa$bounds` data frame.
#'
#' @param sa An object of class \code{sensana}.
#' @param name A character string indicating the name of the bound to be removed.
#'
#' @return An updated \code{sensana} object with the specified bound removed.
#'
#' @examples
#' u <- rnorm(20)
#' xt <- -0.25*u + rnorm(20)
#' xp <- 0.5 * xt + rnorm(20)
#' z <- 0.25 * xt - xp + u + rnorm(20)
#' d <- 0.5 * xt + 0.5 * xp + u + 2*z + rnorm(20)
#' y <- d + 2*xp - xt + z + u + rnorm(20)
#' sa <- sensana(y = y, d = d, x = data.frame(xt = xt, xp = xp), z = z,
#'               dep_x = "xt", indep_x = "xp", alpha = 0.05, quantile = "t")
#'
#' ## Add bound
#' sa <- add_bound(sa, arrow = "UD", kind = "direct", lb = -0.3, ub = 0.3, name = "example_bound")
#' ## Remove bound
#' sa <- remove_bound(sa, name = "example_bound")
#'
#' @export
remove_bound <- function(sa, name) {
  
  if (!inherits(sa, "sensana")) {
    stop("'sa' must be a 'sensana' object.")
  }
  
  if (!is.character(name) || length(name) != 1 ||
      !(name %in% rownames(sa$bounds))) {
    stop("'name' is not a string or not among the names of the specified bounds.")
  }
  
  sa$bounds <- sa$bounds[rownames(sa$bounds) != name, ]
  sa
}



## Comparative bounds for UD and UY
comp_bound_udy <- function(a, xp, xt, z, b, I, J) {
  xpic <- xp[, setdiff(colnames(xp), I)]
  xtpi <- if (is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])
  xpj <- if(dim(xp)[2] < 2) xp else xp[,J]
  
  r2_axj_xxi <- r2(a, xpj, xtpi)
  r2_axic_xxi <- r2(a, xpic, xtpi)
  
  bd <- sqrt(b * r2_axj_xxi / (1 - r2_axic_xxi))
  b_max <- (1 - r2_axic_xxi) / r2_axj_xxi
  c(bd, b_max)
}



## Comparative-d bound for UY
comp_d_bound_uy <- function(y, d, xp, xt, z, b, I, J) {
  xtpi <- if(is.null(I)) cbind(xt, z) else cbind(xt, z, xp[, I])
  xpj <- if(dim(xp)[2] < 2) xp else xp[, J]
  r2_yxj_xxI <- r2(y, xpj, cbind(xtpi, d))
  
  bd <- sqrt(b * r2_yxj_xxI)
  b_max <- 1 / r2_yxj_xxI
  c(bd, b_max)
}


## Comparative bound for ZU
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


## Print warning when the specified b is too large to have an impact on the PIR.
## The bound may still affect the sensitivity interval, however.
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


## Used when bootstrapping: Since comparative bounds are data-dependent,
## they need to be re-computed.
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
