
## Computing and plotting R-contours
#' @export
r_contours <- function(sa, val_interest = 0,
                       comparison_ind = NULL,
                       comparison = c("comp", "comp-d", "informal"),
                       iv_lines = FALSE,
                       bw = FALSE,
                       grid_specs = list(N1 = 100, N2 = 100, N5 = 100),
                       eps = 0.001,
                       n_breaks = 15) {
  
  
  if (!is.numeric(val_interest) || length(val_interest) != 1) {
    stop("'val_interest' must be a number.")
  }
  
  comparison <- match.arg(comparison)
  
  if (!any(bw %in% c(TRUE, FALSE))) {
    stop("'bw' must be a boolean value.")
  }
  
  if (!is.numeric(n_breaks) || length(n_breaks) != 1 ||
      n_breaks <= 1) {
    stop("'n_breaks' must be a at least 2.")
  }
  
  
  plot_list <- r_contours_data(sa, comparison_ind, comparison,
                               iv_lines, grid_specs, eps)
  r_contours_plot(plot_list, val_interest, bw, n_breaks)
}



## Capsule this method

df_plot_val <- function(p1_seq, p2_mat, eval_mat) {
  len_x <- length(p1_seq)
  len_y <- dim(p2_mat)[2]
  
  xtemp <- (p1_seq[2:len_x] + p1_seq[1:(len_x-1)]) / 2
  xmin <- c(2 * p1_seq[1] - xtemp[1], xtemp)
  xmax <- c(xtemp, 2 * p1_seq[len_x] - xtemp[len_x-1])
  ytemp <- (p2_mat[, 2:len_y] + p2_mat[, 1:(len_y-1)]) / 2
  ymin <- cbind(2 * p2_mat[, 1] - ytemp[, 1], ytemp)
  ymax <- cbind(ytemp, 2 * p2_mat[, len_y] - ytemp[, len_y-1])
  
  data.frame(a = rep(p1_seq, len_y),
             b = as.vector(p2_mat),
             val = as.vector(eval_mat),
             xmin = rep(xmin, len_y),
             xmax = rep(xmax, len_y),
             ymin = as.vector(ymin),
             ymax = as.vector(ymax),
             feasible = as.vector(!is.na(p2_mat)))
}


df_comparison_points <- function(sa, comparison_ind, comparison) {
  
  cp_x <- c()
  cp_y <- c()
  cp_labels <- c()
  for (j in names(comparison_ind)) {
    xj <- sa$xp[, j]
    xxjc <- cbind(sa$xt, sa$z, sa$xp[, setdiff(colnames(sa$xp), j)])
    r_dx <- r(sa$d, xj, xxjc)
    for (m in comparison_ind[[j]]) {
      if (comparison == "comp") {
        cd <- sqrt(m) * f(r_dx)
        cy <- sqrt(m) * f(r(sa$y, xj, xxjc))
        r_yd_x <- r(sa$y, sa$d, cbind(sa$xp, sa$xt, sa$z))
        add_cp_x <- cd
        add_cp_y <- (cy - r_yd_x * cd) / sqrt(1 - r_yd_x^2) / sqrt(1 - cd^2)
      } else if (comparison == "comp-d") {
        add_cp_x <- sqrt(m) * f(r_dx)
        r_yx_d <- r(sa$y, xj, cbind(xxjc, sa$d))
        add_cp_y <- (sqrt(m) * sqrt(1 - m * r_dx^2) + sqrt(m) * r_dx^2 / sqrt(1-r_dx^2)) /
          sqrt(1 - m * r_dx^2 / (1 - r_dx^2)) * r_yx_d / sqrt(1 - r_yx_d^2)
      } else { ## comparison == "informal"
        add_cp_x <- sqrt(m) * r_dx
        add_cp_y <- sqrt(m) * r(sa$y, xj, cbind(xxjc, sa$d))
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
  data.frame(cp_x = cp_x, cp_y = cp_y, cp_label = cp_labels)
}


iv_lines_list <- function(sa, p1_seq, p2_mat) {
  c3 <- r(sa$y, sa$z, cbind(sa$xt, sa$xp, sa$d))
  c4 <- r(sa$d, sa$z, cbind(sa$xt, sa$xp))
  c <- - f(c3) / c4
  
  a_min <- p1_seq[1]
  a_max <- p1_seq[length(p1_seq)]
  b_min <- min(p2_mat, na.rm = TRUE)
  b_max <- max(p2_mat, na.rm = TRUE)
  
  fun <- function(a) c / f(a)
  fun_inv <- function(b) sign(b) * c / sqrt(b^2 + c^2)
  b_min_inv <- fun_inv(b_min)
  b_max_inv <- fun_inv(b_max)
  
  list(fun = fun, fun_inv = fun_inv,
       c = c,
       a = c(a_min, a_max),
       b = c(b_min, b_max),
       b_inv = c(b_min_inv, b_max_inv))
}



#' @export
r_contours_data <- function(sa, comparison_ind,
                            comparison = c("comp", "comp-d", "informal"),
                            iv_lines, grid_specs, eps) {
  
  comparison <- match.arg(comparison)
  check_r_contours_data(sa, comparison_ind, iv_lines, grid_specs, eps)
  
  grid_list <- feasible_grid(sa$y, sa$d, sa$xt, sa$xp, sa$z,
                             sa$bounds, grid_specs, TRUE, eps)
  ep <- compute_eval_param(sa$y, sa$d, sa$xt, sa$xp, sa$z)
  
  eval_mat <- eval_on_grid(grid_list$p1_seq, grid_list$p2_mat,
                           ep$beta_ols, ep$sd_y_xzd, ep$sd_d_xz)
  
  df_plot <- df_plot_val(grid_list$p1_seq, grid_list$p2_mat, eval_mat)
  df_cp <- if(is.null(comparison_ind)) NULL else df_comparison_points(sa, comparison_ind, comparison)
  iv_lines_list <- if(iv_lines) iv_lines_list(sa, grid_list$p1_seq, grid_list$p2_mat) else NULL
  
  list(df_plot = df_plot, df_cp = df_cp, ivl = iv_lines_list)
}



#' @importFrom rlang .data
r_contours_plot <- function(plist, val_interest, bw, n_breaks) {

  make_breaks <- function(range, binwidth) {
    signif(pretty(range, n_breaks), 4)
  }

  ## Grid plot
  pl <- ggplot2::ggplot(data = plist$df_plot[plist$df_plot$feasible,]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data$val, colour = .data$val),
                       alpha = 1, na.rm = TRUE)
  
  if (bw) {
    pl <- pl + ggplot2::scale_fill_steps(low = "#5A5A5A", high = "#F5F5F5",
                                         breaks = make_breaks,
                                         aesthetics = c("fill", "colour"))
  } else {
    pl <- pl + ggplot2::scale_fill_steps2(midpoint = val_interest,
                                          breaks = make_breaks,
                                          aesthetics = c("fill", "colour"))
  }
  
  ## IV-lines
  if (!is.null(iv_lines_list)) {
    if (plist$ivl$c > 0) {
      if (plist$ivl$b[2] > 0 && plist$ivl$b_inv[2] < plist$ivl$a[2]) {
        pl <- pl + ggplot2::stat_function(fun = plist$ivl$fun,
                                          xlim = c(plist$ivl$b_inv[2], plist$ivl$a[2]),
                                          size = 1)
      }
      if (plist$ivl$b[1] < 0 && plist$ivl$b_inv[1] > plist$ivl$a[1]) {
        pl <- pl + ggplot2::stat_function(fun = plist$ivl$fun,
                                          xlim = c(plist$ivl$a[1], plist$ivl$b_inv[1]),
                                          size = 1)
      }
    } else {
      if (plist$ivl$b[2] > 0 && plist$ivl$b_inv[2] > plist$ivl$a[1]) {
        pl <- pl + ggplot2::stat_function(fun = plist$ivl$fun,
                                          xlim = c(plist$ivl$a[1], plist$ivl$b_inv[2]),
                                          size = 1)
      }
      if (plist$ivl$b[1] < 0 && plist$ivl$b_inv[1] < plist$ivl$a[2]) {
        pl <- pl + ggplot2::stat_function(fun = plist$ivl$fun,
                                          xlim = c(plist$ivl$b_inv[1], plist$ivl$a[2]),
                                          size = 1)
      }
    }
  }

  ## Comparison points
  if (!is.null(plist$df_cp) && dim(plist$df_cp)[1] > 0) {
    pl <- pl +
      ggplot2::geom_point(data = plist$df_cp,
                 mapping = ggplot2::aes(x = .data$cp_x, y = .data$cp_y),
                 col = "black") +
      ggrepel::geom_text_repel(data = plist$df_cp,
                               mapping = ggplot2::aes(x = .data$cp_x,
                                                      y = .data$cp_y,
                                                      label = .data$cp_label),
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
