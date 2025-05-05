check_sensana <- function(y, d, indep_x, dep_x, x, z, quantile, alpha) {
  
  if (!is.vector(y) || !is.numeric(y)) {
    stop("'y' must be a numeric vector.")
  }
  
  if (!is.vector(d) || !is.numeric(d)) {
    stop("'d' must be a numeric vector.")
  }
  
  if (!is.null(indep_x) && !is.character(indep_x)) {
    stop("'indep_x must be NULL or a character vector.")
  }
  
  if (!is.null(dep_x) && !is.character(dep_x)) {
    stop("'dep_x must be NULL or a character vector.")
  }
  
  if (!is.null(x) && (!is.data.frame(x) || !is.numeric(as.matrix(x)))) {
    stop("'x' must be NULL or a numeric data frame.")
  }
  
  if (!is.null(z) && (!is.vector(z) || !is.numeric(z))) {
    stop("'z' must be NULL or a numeric vector.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >=1) {
    stop("'alpha' must be a real-number within (0,1).")
  }
  
  if (any(indep_x %in% dep_x) || any(dep_x %in% indep_x)) {
    stop("'indep_x' and 'dep_x' must not overlap.")
  }
  
  if (!all(colnames(x) %in% c(indep_x, dep_x)) ||
      !all(c(indep_x, dep_x) %in% colnames(x))) {
    stop("'indep_x' and 'dep_x' must jointly contain all column names of x.")
  }
  
  if (length(y) != length(d) ||
      (!is.null(z) && length(z) != length(d)) ||
      (!is.null(x) && dim(x)[1] != length(d))) {
    stop("The length/dimension of 'y', 'd', 'x' and 'z' must align.")
  }
  
}


check_add_bound <- function(sa, arrow, kind, lb, ub, b, I, J,
                            name, print_warning) {
  
  if (!inherits(sa, "sensana")) {
    stop("'sa' must be a 'sensana' object.")
  }
  
  if (kind == "direct" && (!is.numeric(lb) || !is.numeric(ub) ||
                           length(lb) != 1 || length(ub) != 1 ||
                           lb < -1 || ub > 1 || lb > ub)) {
    stop("'lb' and 'ub' must form an interval within [-1,1].")
  }
  
  if (kind != "direct" && (!is.numeric(b) || length(b) != 1 || b <= 0)) {
    stop("'b' must be a positive number.")
  }
  
  if (!is.character(name)) {
    stop("'name' must be a string.")
  }
  
  if (!any(print_warning %in% c(TRUE, FALSE))) {
    stop("'print_warning' must be a boolean value.")
  }
  
  if (kind == "comparative-d" && arrow != "UY") {
    stop("'comparative-d' bounds are only possible for the arrow 'UY'.")
  }
  
  if (arrow %in% c("UD", "UY") && kind != "direct") {
    if (!is.null(I) && !is.character(I)) {
      stop("'I' must be NULL, a string or a vector of strings.")
    }
    
    if (!is.character(J)) {
      stop("'J' must be a string or a vector of strings.")
    }
    
    if (any(J %in% I) || any(I %in% J)) {
      stop("'I' and 'J' must not overlap.")
    }
    
    if (!all(c(J,I) %in% colnames(sa$xp))) {
      stop("'I' and 'J' must be subsets of the conditionall independent covariates.")
    }
    
  }
  if (arrow %in% c("ZU", "ZY") && kind != "direct") {
    if (!is.character(J) || length(J) != 1) {
      stop("'J' must be one string for bounds on 'ZU' and 'ZY'.")
    }
    if (!(J %in% colnames(sa$xp))) {
      stop("'J' must be the name of a conditionally independent covariate.")
    }
  }
  
}


check_pir <- function(sa, grid_specs, eps) {
  if (!inherits(sa, "sensana")) {
    stop("'sa' must be a 'sensana' object.")
  }
  
  if (!is.list(grid_specs) ||
      !all(names(grid_specs) == c("N1", "N2", "N5")) ||
      !is.numeric(grid_specs$N1) || !is.numeric(grid_specs$N2) ||
      !is.numeric(grid_specs$N5) || length(grid_specs$N1) != 1 ||
      length(grid_specs$N2) != 1 || length(grid_specs$N5) != 1) {
    stop("'grid_specs' must be list of 3 numbers: N1, N2, N5")
  }
  
  if (any(unlist(grid_specs) < 2)) {
    stop("The values in 'grid_specs' must be at least 2.")
  }
  
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0 || eps > 0.1) {
    stop("'eps' must be a real-number within (0,0.1].")
  }
}


check_sensint <- function(sa, alpha, boot_procedure, boot_samples,
                          grid_specs, eps, parallel, ncpus) {
  check_pir(sa, grid_specs, eps)
  
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >=1) {
    stop("'alpha' must be a real-number within (0,1).")
  }
  
  if (!is.numeric(boot_samples) || length(boot_samples) != 1 || boot_samples < 2) {
    stop("'boot_samples' must be a real-number strictly larger than 1.")
  }
  
  if (!is.character(boot_procedure) ||
      !all(boot_procedure %in% c("norm", "basic", "stud", "perc", "bca"))) {
    stop(paste0("'boot_procedure' must be a character or vector of characters",
                " taking the values: 'norm', 'basic', 'stud', 'perc', 'bca'."))
  }
}


check_b_contours_data <- function(sa, pir_lower, bound1, range1, bound2, range2,
                                  grid_specs_b, grid_specs, eps) {
  
  if (!inherits(sa, "sensana")) {
    stop("'sa' must be a 'sensana' object.")
  }
  
  if (!any(pir_lower %in% c(TRUE, FALSE))) {
    stop("'pir_lower' must be a boolean value.")
  }
  
  if (!is.character(bound1) || length(bound1) != 1 ||
      !(bound1 %in% rownames(sa$bounds))) {
    stop("'bound1' is not a string or not among the names of the specified bounds.")
  }
  
  if (!is.character(bound2) || length(bound2) != 1 ||
      !(bound2 %in% rownames(sa$bounds))) {
    stop("'bound2' is not a string or not among the names of the specified bounds.")
  }
  
  if (!is.numeric(range1) || length(range1) != 2 ||
      range1[1] <= 0 || range1[1] >= range1[2]) {
    stop("'range1' must be a vector of two ordered positive numbers.")
  }
  
  if (!is.numeric(range2) || length(range2) != 2 ||
      range2[1] <= 0 || range2[1] >= range2[2]) {
    stop("'range2' must be a vector of two ordered positive numbers.")
  }
  
  if (!is.list(grid_specs_b) ||
      !all(names(grid_specs_b) == c("N1b", "N2b")) ||
      !is.numeric(grid_specs_b$N1b) || !is.numeric(grid_specs_b$N2b) ||
      length(grid_specs_b$N1b) != 1 || length(grid_specs_b$N2b) != 1) {
    stop("'grid_specs_b' must be list of 2 numbers: N1b, N2b")
  }
  
  if (any(unlist(grid_specs_b) < 2)) {
    stop("The values in 'grid_specs_b' must be at least 2.")
  }
  
  if (!is.list(grid_specs) ||
      !all(names(grid_specs) == c("N1", "N2", "N5")) ||
      !is.numeric(grid_specs$N1) || !is.numeric(grid_specs$N2) ||
      !is.numeric(grid_specs$N5) || length(grid_specs$N1) != 1 ||
      length(grid_specs$N2) != 1 || length(grid_specs$N5) != 1) {
    stop("'grid_specs' must be list of 3 numbers: N1, N2, N5")
  }
  
  if (any(unlist(grid_specs) < 2)) {
    stop("The values in 'grid_specs' must be at least 2.")
  }
  
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0 || eps > 0.1) {
    stop("'eps' must be a real-number within (0,0.1].")
  }
}


check_r_contours_data <- function(sa, comparison_ind,
                                  iv_lines, grid_specs, eps) {
  
  if (!inherits(sa, "sensana")) {
    stop("'sa' must be a 'sensana' object.")
  }
  
  if (!is.null(comparison_ind) &&
      (!is.list(comparison_ind) ||
       !all(names(comparison_ind) %in% colnames(sa$xp)) ||
       !all(is.numeric(unlist(comparison_ind))) ||
       !all(unlist(comparison_ind) > 0))) {
    stop(paste0("'comparison_ind' must be NULL or a list, where each element ",
                "has the name of a conditionally independent covariate ",
                "and contains a (vector of) positive number(s)."))
  }
  
  if (!any(iv_lines %in% c(TRUE, FALSE))) {
    stop("'iv_lines' must be a boolean value.")
  }
  
  if (!is.list(grid_specs) ||
      !all(names(grid_specs) == c("N1", "N2", "N5")) ||
      !is.numeric(grid_specs$N1) || !is.numeric(grid_specs$N2) ||
      !is.numeric(grid_specs$N5) || length(grid_specs$N1) != 1 ||
      length(grid_specs$N2) != 1 || length(grid_specs$N5) != 1) {
    stop("'grid_specs' must be list of 3 numbers: N1, N2, N5")
  }
  
  if (any(unlist(grid_specs) < 2)) {
    stop("The values in 'grid_specs' must be at least 2.")
  }
  
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0 || eps > 0.1) {
    stop("'eps' must be a real-number within (0,0.1].")
  }
}