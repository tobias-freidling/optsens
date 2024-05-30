
## Computing and plotting R-contours
#' @export
r_contours <- function(sa, val_interest = 0,
                       comparison_ind = NULL,
                       mult_r2 = TRUE,
                       comparison = "comp", ## other values "comp-d" and "informal"
                       iv_lines = FALSE,
                       grid_specs = list(num_x = 100,
                                         num_y = 100,
                                         num_z = 100),
                       print_warning = FALSE,
                       eps = 0.001) {
  plot_list <- r_contours_data(sa, comparison_ind,
                               mult_r2, comparison,
                               iv_lines, grid_specs,
                               print_warning, eps)
  r_contours_plot(plot_list, val_interest)
}



#' @export
r_contours_data <- function(sa, comparison_ind, mult_r2, comparison,
                            iv_lines, grid_specs, print_warning, eps) {
  list2env(sa, environment())

  grid_list <- feasible_grid(y, d, xt, xp, z, bounds, grid_specs, TRUE,
                             print_warning, eps)
  list2env(grid_list, environment())
  eval_param <- compute_eval_param(y, d, xt, xp, z)
  eval_mat <- eval_on_grid(a_seq, b_mat,
                           eval_param[1], eval_param[2], eval_param[3])

  ## Plot info
  len_x <- length(a_seq)
  len_y <- dim(b_mat)[2]

  xtemp <- (a_seq[2:len_x] + a_seq[1:(len_x-1)]) / 2
  xmin <- c(2 * a_seq[1] - xtemp[1], xtemp)
  xmax <- c(xtemp, 2 * a_seq[len_x] - xtemp[len_x-1])
  ytemp <- (b_mat[, 2:len_y] + b_mat[, 1:(len_y-1)]) / 2
  ymin <- cbind(2 * b_mat[, 1] - ytemp[, 1], ytemp)
  ymax <- cbind(ytemp, 2 * b_mat[, len_y] - ytemp[, len_y-1])

  df_plot <- data.frame(a = rep(a_seq, len_y),
                        b = as.vector(b_mat),
                        val = as.vector(eval_mat),
                        xmin = rep(xmin, len_y),
                        xmax = rep(xmax, len_y),
                        ymin = as.vector(ymin),
                        ymax = as.vector(ymax),
                        feasible = as.vector(!is.na(b_mat)))

  ## Comparison points
  cp_x <- c()
  cp_y <- c()
  cp_labels <- c()
  for (j in names(comparison_ind)) {
    xj <- xp[, j]
    xxjc <- cbind(xt, z, xp[, setdiff(colnames(xp), j)])
    r_dx <- r(d, xj, xxjc)

    for (m in comparison_ind[[j]]) {
      fac <- if(mult_r2) sqrt(m) else m

      if (comparison == "comp") {
        r_yx <- r(y, xj, xxjc)
        r_yd_x <- r(y, d, cbind(xp, xt, z))
        cd <- fac * r_dx / sqrt(1 - r_dx^2)
        cy <- fac * r_yx / sqrt(1 - r_yx^2)
        add_cp_x <- cd
        add_cp_y <- (cy - r_yd_x * cd) /
          sqrt(1 - r_yd_x^2) / sqrt(1 - cd^2)
      } else if (comparison == "comp-alt") {
        r_yx_d <- r(y, xj, cbind(xxjc, d))
        add_cp_x <- fac * r_dx / sqrt(1 - r_dx^2)
        add_cp_y <- fac / sqrt(1 - (1 + fac^2) * r_dx^2) * r_yx_d / sqrt(1 - r_yx_d^2)
      } else if (comparison == "comp-d") {
        add_cp_x <- fac * r_dx / sqrt(1 - r_dx^2)
        r_yx_d <- r(y, xj, cbind(xxjc, d))
        add_cp_y <- (fac * sqrt(1 - fac^2 * r_dx^2) + fac * r_dx^2 / sqrt(1-r_dx^2)) /
          sqrt(1 - fac^2 * r_dx^2 / (1 - r_dx^2)) * r_yx_d / sqrt(1 - r_yx_d^2)
      } else {
        add_cp_x <- fac * r_dx
        add_cp_y <- fac * r(y, xj, cbind(xxjc, d))
      }

      if (abs(add_cp_x) >= 1 || abs(add_cp_y) >= 1) {
        stop(paste0("The multiplier ", m, " is too large ",
                    "for the covariate ", j, "."))
      }

      add_label <- if(m == 1) j else paste0(m, "x ", j)
      cp_x <- c(cp_x, add_cp_x)
      cp_y <- c(cp_y, add_cp_y)
      cp_labels <- c(cp_labels, add_label)
    }
  }
  df_cp <- data.frame(cp_x = cp_x, cp_y = cp_y,
                      cp_label = cp_labels)

  ## IV-lines
  if (iv_lines) {
    c3 <- r(y, z, cbind(xt, xp, d))
    c4 <- r(d, z, cbind(xt, xp))
    c <- -c3 / sqrt(1-c3^2) / c4

    a_min <- a_seq[1]
    a_max <- a_seq[length(a_seq)]
    b_min <- min(b_mat, na.rm = TRUE)
    b_max <- max(b_mat, na.rm = TRUE)

    fun <- function(a) c * sqrt(1-a^2) / a
    fun_inv <- function(b) sign(b) * c / sqrt(b^2 + c^2)
    b_min_inv <- fun_inv(b_min)
    b_max_inv <- fun_inv(b_max)

    iv_lines_list <- list(fun = fun, fun_inv = fun_inv,
                          c = c,
                          a = c(a_min, a_max),
                          b = c(b_min, b_max),
                          b_inv = c(b_min_inv, b_max_inv))
  } else {
    iv_lines_list <- NULL
  }

  list(df_plot = df_plot, df_cp = df_cp,
       iv_lines_list = iv_lines_list)
}



#' @export
r_contours_plot <- function(plot_list, val_interest) {

  list2env(plot_list, environment())

  make_breaks <- function(range, binwidth) {
    signif(pretty(range, 15), 4)
  }

  ## Grid plot
  pl <- ggplot2::ggplot(subset(df_plot, feasible)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax,
                                    fill = val, colour = val),
                       alpha = 1, na.rm = TRUE) +
    ggplot2::scale_fill_steps2(midpoint = val_interest,
                               breaks = make_breaks,
                               aesthetics = c("fill", "colour"))

  ## IV-lines
  if (!is.null(iv_lines_list)) {
    list2env(plot_list$iv_lines_list, environment())
    if (c > 0) {
      if (b[2] > 0 && b_inv[2] < a[2]) {
        pl <- pl + ggplot2::stat_function(fun = fun,
                                          xlim = c(b_inv[2], a[2]),
                                          size = 1)
      }
      if (b[1] < 0 && b_inv[1] > a[1]) {
        pl <- pl + ggplot2::stat_function(fun = fun,
                                          xlim = c(a[1], b_inv[1]),
                                          size = 1)
      }
    } else {
      if (b[2] > 0 && b_inv[2] > a[1]) {
        pl <- pl + ggplot2::stat_function(fun = fun,
                                          xlim = c(a[1], b_inv[2]),
                                          size = 1)
      }
      if (b[1] < 0 && b_inv[1] < a[2]) {
        pl <- pl + ggplot2::stat_function(fun = fun,
                                          xlim = c(b_inv[1], a[2]),
                                          size = 1)
      }
    }

  }

  ## Comparison points
  if (dim(df_cp)[1] > 0) {
    list2env(plot_list$df_cp, environment())
    pl <- pl +
      ggplot2::geom_point(data = df_cp,
                 mapping = ggplot2::aes(x = cp_x, y = cp_y),
                 col = "black") +
      ggrepel::geom_text_repel(data = df_cp,
                               mapping = ggplot2::aes(x = cp_x, y = cp_y,
                                             label = cp_label),
                               col = "black", size = 3.5)
  }

  xlab <- expression("R"["D~U|X,Z"])
  ylab <- expression("R"["Y~U|X,Z,D"])
  title <- "R-sensitivity contours"

  pl <- pl +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text = ggplot2::element_text(size = 11),
                   axis.title = ggplot2::element_text(size = 13),
                   title = ggplot2::element_text(size = 13),
                   legend.key.size = grid::unit(1.5, 'cm'),
                   legend.title = ggplot2::element_blank())
  pl
}
