
## Computing and plotting b-contours
#' @export
b_contours <- function(sa, pir_lower, bound1, range1, bound2, range2,
                       val_interest = 0, grid_specs_b = list(N1b = 20, N2b =20),
                       grid_specs = list(N1 = 100, N2 = 100, N5 = 100),
                       eps = 0.001, bw = FALSE) {
  
  if (!is.numeric(val_interest) || length(val_interest) != 1) {
    stop("'val_interest' must be a number.")
  }
  
  if (!any(bw %in% c(TRUE, FALSE))) {
    stop("'bw' must be a boolean value.")
  }

  data <- b_contours_data(sa, pir_lower, bound1, range1, bound2, range2,
                          grid_specs_b, grid_specs, eps)

  xlab <- sa$bounds[bound1, "arrow"]
  ylab <- sa$bounds[bound2, "arrow"]
  b1_sa <- sa$bounds[bound1, "b"]
  b2_sa <- sa$bounds[bound2, "b"]

  b_contours_plot(data, val_interest, pir_lower, b1_sa, b2_sa,
                  range1, range2, xlab, ylab, bw)
}



#' @export
b_contours_data <- function(sa, pir_lower, bound1, range1, bound2, range2,
                            grid_specs_b, grid_specs, eps){
  
  check_b_contours_data(sa, pir_lower, bound1, range1, bound2, range2,
                        grid_specs_b, grid_specs, eps)

  b1_seq <- seq(range1[1], range1[2], length = grid_specs_b$N1b)
  b2_seq <- seq(range2[1], range2[2], length = grid_specs_b$N2b)

  data_mat <- cbind(sa$y, sa$d, sa$z, sa$xp, sa$xt)
  if (is.null(sa$z)) {
    colnames(data_mat)[c(1,2)] <- c("y", "d")
  } else {
    colnames(data_mat)[c(1,2,3)] <- c("y", "d", "z")
  }
  indices <- 1:length(sa$y)
  indep_x <- if (is.null(sa$xp)) NULL else colnames(sa$xp)
  dep_x <- if (is.null(sa$xt)) NULL else colnames(sa$xt)

  pir_b1b2 <- function(b1, b2) {
    bounds_b1b2 <- sa$bounds
    bounds_b1b2[bound1, "b"] <- b1
    bounds_b1b2[bound2, "b"] <- b2
    
    pir <- one_opt(data_mat, indices, bounds_b1b2, grid_specs, indep_x, dep_x, eps)
    if (pir_lower) pir[1] else pir[2]
  }

  pir_val <- outer(b1_seq, b2_seq, Vectorize(pir_b1b2))

  cbind(expand.grid(x = b1_seq, y = b2_seq), z = as.vector(pir_val))
}



#' @importFrom rlang .data
b_contours_plot <- function(data, val_interest, pir_lower,
                            b1_sa, b2_sa, range1, range2,
                            xlab, ylab, bw) {

  title <- "b-sensitivity contours: "
  title <- if (pir_lower) paste0(title, "lower point PIR") else paste0(title, "upper point PIR")

  text_point <- paste0("(", b1_sa, ", ", b2_sa, ")")

  make_breaks <- function(range, binwidth) {
    signif(pretty(range, 15), 4)
  }
  make_breaks_ex <- function(range, binwidth) {
    b <- make_breaks(range, binwidth)
    b[b != val_interest]
  }

  pl <- ggplot2::ggplot(data, ggplot2::aes(x = .data$x, y = .data$y, na.rm = TRUE)) +
    metR::geom_contour_fill(ggplot2::aes(z = .data$z,
                                         fill = ggplot2::after_stat(.data$level)),
                            breaks = make_breaks,
                            show.legend = FALSE) +
    metR::geom_contour2(ggplot2::aes(z = .data$z,
                            label = ggplot2::after_stat(.data$level)),
                        breaks = make_breaks_ex,
                        col = "black",
                        label_size = 3,
                        size = 0.1) +
    metR::geom_contour2(ggplot2::aes(z = .data$z,
                                     label = ggplot2::after_stat(.data$level)),
                        breaks = val_interest,
                        linewidth = 0.8,
                        col = "black")
  
  if (bw) {
    pl <- pl + metR::scale_fill_discretised(low = "#5A5A5A", high = "#F5F5F5")
  } else {
    pl <- pl + metR::scale_fill_divergent_discretised(midpoint = val_interest)
  }
  
  pl <- pl +
    ggplot2::geom_point(data = data.frame(x = b1_sa,
                                          y = b2_sa),
                        mapping = ggplot2::aes(x = .data$x, y = .data$y),
                        col = "black") +
    ggrepel::geom_text_repel(data = data.frame(x = b1_sa,
                                               y = b2_sa,
                                               label = text_point),
                             mapping = ggplot2::aes(x = .data$x,
                                                    y = .data$y,
                                                    label = .data$label),
                             col = "black",
                             size = 3.5,
                             xlim = range1,
                             ylim = range2) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text = ggplot2::element_text(size = 11),
                   axis.title = ggplot2::element_text(size = 13),
                   title = ggplot2::element_text(size = 12))
  pl
}
