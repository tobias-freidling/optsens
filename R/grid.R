
## Computing the constants used in the grid search algorithm
constants <- function(y, d, xt, xp, z, bounds) {
  xtz <- cbind(xt, z)
  x <- cbind(xt, xp)
  
  ## OLS-related
  c1 <- r(y, d, cbind(xtz, xp))
  c2 <- c3 <- c4 <- lb_p4 <- ub_p4 <- c()
  for (i in seq_len(nrow(bounds))) {
    if (bounds[i, "arrow"] == "UY" &&
        bounds[i, "kind"] == "comparative-d") {
      I <- unlist(bounds[i, "I"])
      xpic <- xp[, setdiff(colnames(xp), I)]
      xxi <- if (is.null(I)) xtz else cbind(xtz, xp[, I])
      c2 <- c(c2, r(y, d, xxi))
      c3 <- c(c3, r(d, xpic, xxi))
      c4 <- c(c4, r(y, xpic, xxi))
      lb_p4 <- c(lb_p4, bounds[i, "lb"])
      ub_p4 <- c(ub_p4, bounds[i, "ub"])
    }
  }
  
  if (is.null(c2)) {
    c2 <- c3 <- c4 <- lb_p4 <- ub_p4 <- NA 
  }
  
  ## TSLS-related
  if (is.null(z)) {
    c5 <- c6 <- NA
  } else {
    c5 <- r(d, z, x)
    c6 <- r(y, z, cbind(x, d))
  }
  
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
  
  if (is.null(c7)) {
    c7 <- b7 <- NA
  }
  
  list(c1 = c1, c2 = c2, c3 = c3, c4 = c4, c5 = c5, c6 = c6, c7 = c7, b7 = b7,
       lb_p4 = lb_p4, ub_p4 = ub_p4)
}




## Extracting and gathering the upper and lower bounds from the bounds data frame
static_bounds <- function(bounds, eps) {
  
  ## Bounds on p1 (short for psi1)
  bp1 <- bounds[bounds$arrow == "UD", c("lb", "ub")]
  lb_p1s <- max(-1 + eps, bp1$lb)
  ub_p1s <- min(1 - eps, bp1$ub)
  
  ## Bounds on p2
  bp2 <- bounds[(bounds$arrow == "UY") & (bounds$kind == "direct"), c("lb", "ub")]
  lb_p2s <- max(-1 + eps, bp2$lb)
  ub_p2s <- min(1 - eps, bp2$ub)
  
  ## Bounds on p3
  bp3 <- bounds[(bounds$arrow == "UY") & (bounds$kind == "comparative"), c("lb", "ub")]
  lb_p3s <- max(-1 + eps, bp3$lb)
  ub_p3s <- min(1 - eps, bp3$ub)
  
  ## for TSLS constraints
  if (any(bounds[,"arrow"] %in% c("ZU", "ZY"))) {
    ## Bound on p6
    bp6 <- bounds[bounds$arrow == "ZU", c("lb", "ub")]
    lb_p6s <- max(-1 + eps, bp6$lb)
    ub_p6s <- min(1 - eps, bp6$ub)
    
    ## Bounds on p7
    bp7 <- bounds[(bounds$arrow == "ZY") & (bounds$kind == "direct"), c("lb", "ub")]
    lb_p7s <- max(-1 + eps, bp7$lb)
    ub_p7s <- min(1 - eps, bp7$ub)
  } else {
    lb_p6s <- ub_p6s <- lb_p7s <- ub_p7s <- NA
  }
  
  list(lb_p1s = lb_p1s, ub_p1s = ub_p1s, lb_p2s = lb_p2s, ub_p2s = ub_p2s,
       lb_p3s = lb_p3s, ub_p3s = ub_p3s, lb_p6s = lb_p6s, ub_p6s = ub_p6s,
       lb_p7s = lb_p7s, ub_p7s = ub_p7s)
}



#' Compute a grid of feasible points.
#'
#' Constructs a grid of \eqn{(\psi_1,\psi_2)} points that are (approximately) contained in the
#' sensitivity model \eqn{\Psi} which is defined via the user-specified bounds.
#' This method should only be used to directly access the grid-search algorithm.
#' For inference, see [pir()] and [sensint()]; for visualizations, see [b_contours()]
#' and [r_contours()].
#'
#'
#' @param y A numeric vector of _centered_ outcomes.
#' @param d A numeric vector of _centered_ treatments.
#' @param xt A matrix of _centered_ not conditionally independent covariates.
#' @param xp A matrix of _centered_ conditionally independent covariates.
#' @param z A numeric vector of _centered_ instruments.
#' @param bounds A data frame containing bounds according to the format used by
#'   the class `"sensana"`.
#' @param grid_specs A named list of the three numeric values `N1`, `N2` and `N5` specifying
#'   the number of points considered in the grid search for each of the three dimensions.
#' @param full_grid Logical; if \code{TRUE}, provides the full grid of
#'   \eqn{(\psi_1,\psi_2)} values. If \code{FALSE}, provides only the smallest and
#'   largest \eqn{\psi_2} value for each \eqn{\psi_1}. (This suffices for estimating
#'   the partially identified range.)
#' @param eps A small numeric number to bound sensitivity parameters away from `1`
#'   and `-1`. Default is \code{0.001}.
#'
#' @return A list containing two elements:
#' \describe{
#' \item{p1_seq  }{A numeric vector of length \code{N1} containing the grid points for
#' the sensitivity parameter \eqn{\psi_1}.}
#' \item{p2_mat  }{A numeric matrix with dimension \code{(N1, N2)} if \code{full_grid == TRUE} and
#' \code{(N1, 2)} otherwise. Row \code{i} of the matrix contains the grid points for \eqn{\psi_2}
#' belonging to the \code{i}-th grid point of \eqn{\psi_1}, i.e., the \code{i}-th entry in \code{p1_seq}.}
#' }
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
#' grid <- feasible_grid(y = sa$y, d = sa$d, xt = sa$xt, xp = sa$xp, z = sa$z,
#'                       bounds = sa$bounds,
#'                       grid_specs = list(N1 = 100, N2 = 100, N5 = 100),
#'                       full_grid = TRUE, eps = 0.001)
#'
#' @export
feasible_grid <- function(y, d, xt, xp, z, bounds, grid_specs,
                          full_grid, eps = 0.001) {
  
  sb <- static_bounds(bounds, eps)
  co <- constants(y, d, xt, xp, z, bounds)
  
  exist_tsls_bound <- any(bounds[,"arrow"] %in% c("ZU", "ZY"))
  exist_comp_uy_bound <- dim(bounds[(bounds$arrow == "UY") &
                                      (bounds$kind != "direct"),])[1] > 0
  
  ## Calling the Rcpp function
  res <- grid_search(grid_specs$N1, grid_specs$N2, grid_specs$N5,
                     full_grid,
                     sb$lb_p1s, sb$ub_p1s,
                     sb$lb_p2s, sb$ub_p2s,
                     sb$lb_p3s, sb$ub_p3s,
                     sb$lb_p6s, sb$ub_p6s,
                     sb$lb_p7s, sb$ub_p7s,
                     co$lb_p4, co$ub_p4,
                     exist_comp_uy_bound,
                     exist_tsls_bound,
                     co$c1, co$c2, co$c3, co$c4,
                     co$c5, co$c6, co$c7, co$b7)
  res
}





## Evaluating objective, i.e. causal effect, over the grid
eval_on_grid <- function(p1_seq, p2_mat, beta_ols, sd_y_xzd, sd_d_xz) {
  p1_mat <- matrix(rep(p1_seq, dim(p2_mat)[2]), length(p1_seq), dim(p2_mat)[2])
  beta_ols - p2_mat * f(p1_mat) * sd_y_xzd / sd_d_xz
}
