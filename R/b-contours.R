
## Computing and plotting b-contours
#' @export
b_contours <- function(sa, pir_lower, bound1, range1, bound2, range2,
                       val_interest = 0, grid_specs_b = list(num1 = 20, num2 =20),
                       grid_specs = list(num_x = 100, num_y = 100, num_z = 100),
                       print_warning = FALSE, eps = 0.001) {

  data <- b_contours_data(sa, pir_lower, bound1, range1, bound2, range2,
                          grid_specs_b, grid_specs, print_warning, eps)

  xlab <- sa$bounds[bound1, "arrow"]
  ylab <- sa$bounds[bound2, "arrow"]
  b1_sa <- sa$bounds[bound1, "b"]
  b2_sa <- sa$bounds[bound2, "b"]

  b_contours_plot(data, val_interest, pir_lower, b1_sa, b2_sa,
                  range1, range2, xlab, ylab)
}



#' @export
b_contours_data <- function(sa, pir_lower, bound1, range1, bound2, range2,
                            grid_specs_b, grid_specs, print_warning, eps){
  list2env(sa, environment())

  b1_seq <- seq(range1[1], range1[2], length = grid_specs_b$num1)
  b2_seq <- seq(range2[1], range2[2], length = grid_specs_b$num2)

  data_mat <- cbind(y, d, z, xp, xt)
  if (is.null(z)) {
    colnames(data_mat)[c(1,2)] <- c("y", "d")
  } else {
    colnames(data_mat)[c(1,2,3)] <- c("y", "d", "z")
  }
  indices <- 1:length(y)
  indep_x <- if (is.null(xp)) NULL else colnames(xp)
  dep_x <- if (is.null(xt)) NULL else colnames(xt)

  pir_b1b2 <- function(b1, b2) {
    bounds_b1b2 <- bounds
    bounds_b1b2[bound1, "b"] <- b1
    bounds_b1b2[bound2, "b"] <- b2
    pir <- one_opt(data_mat, indices, bounds_b1b2, grid_specs, indep_x, dep_x,
                   print_warning, eps, NULL)
    if (pir_lower) pir[1] else pir[2]
  }

  pir_val <- outer(b1_seq, b2_seq, Vectorize(pir_b1b2))

  cbind(expand.grid(x = b1_seq, y = b2_seq), z = as.vector(pir_val))
}


#' @export
b_contours_plot <- function(data, val_interest, pir_lower,
                            b1_sa, b2_sa, range1, range2,
                            xlab, ylab) {

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

  pl <- ggplot2::ggplot(data, ggplot2::aes(x, y)) +
    metR::geom_contour_fill(ggplot2::aes(z = z, fill = stat(level)),
                            breaks = make_breaks,
                            show.legend = FALSE,
                            na.rm = TRUE) +
    metR::geom_contour2(ggplot2::aes(z = z,
                            label = stat(level)),
                        breaks = make_breaks_ex,
                        col = "black",
                        label_size = 3,
                        size = 0.1,
                        na.rm = TRUE) +
    metR::geom_contour2(ggplot2::aes(z = z, label = stat(level)),
                        breaks = val_interest,
                        size = 0.8,
                        col = "black",
                        na.rm = TRUE) +
    metR::scale_fill_divergent_discretised(midpoint = val_interest) +
    ggplot2::geom_point(data = data.frame(x = b1_sa,
                                          y = b2_sa),
                        mapping = ggplot2::aes(x, y), col = "black") +
    ggrepel::geom_text_repel(data = data.frame(x = b1_sa,
                                               y = b2_sa,
                                               label = text_point),
                             mapping = ggplot2::aes(x, y, label = label),
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
