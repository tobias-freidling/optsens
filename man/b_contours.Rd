% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/b-contours.R
\name{b_contours}
\alias{b_contours}
\title{Generate a B-Contour Plot}
\usage{
b_contours(
  sa,
  pir_lower,
  bound1,
  range1,
  bound2,
  range2,
  val_interest = 0,
  grid_specs_b = list(N1b = 20, N2b = 20),
  grid_specs = list(N1 = 100, N2 = 100, N5 = 100),
  eps = 0.001,
  bw = FALSE
)
}
\arguments{
\item{sa}{An object of the class \code{sensana}.}

\item{pir_lower}{Logical; if \code{TRUE}, the lower end of the PIR is plotted.
Otherwise, the upper end is plotted.}

\item{bound1}{Character. The name of the first comparative bound.}

\item{range1}{Numeric vector of length 2 specifying the lower and upper ends
of the range of the first b-factor.}

\item{bound2}{Character. The name of the second comparative bound.}

\item{range2}{Numeric vector of length 2 specifying the lower and upper ends
of the range of the second b-factor.}

\item{val_interest}{Numeric; the value of interest for the user; Default is \code{0}.}

\item{grid_specs_b}{A named list of the two numeric values \code{N1b} and \code{N2b},
specifying the resolution of the contour plot, i.e. the number of b-factors
for which the PIR is computed.}

\item{grid_specs}{A named list of the three numeric values \code{N1}, \code{N2} and \code{N5} specifying
the number of points considered in the grid search for each of the three dimensions.}

\item{eps}{A small numeric number to bound sensitivity parameters away from \code{1}
and \code{-1}. Default is \code{0.001}.}

\item{bw}{Logical; whether to generate a black-and-white plot. Default is \code{FALSE}.}
}
\value{
An object of the class \code{ggplot} that can be printed.
}
\description{
This method can be used for sensitivity models that contain (at least) two
comparative bounds and visualizes how the b-factors in the comparative bounds
affect the estimated partially identified range (PIR). Users specify a range of values
for the two b-factors and the method plots the contour lines for the upper or lower
end of the PIR.
}
\examples{
set.seed(123)
u <- rnorm(20)
xt <- -0.25*u + rnorm(20)
xp <- 0.5 * xt + rnorm(20)
z <- 0.25 * xt - xp + u + rnorm(20)
d <- 0.5 * xt + 0.5 * xp + u + 2*z + rnorm(20)
y <- d + 2*xp - xt + z + u + rnorm(20)
sa <- sensana(y = y, d = d, x = data.frame(xt = xt, xp = xp), z = z,
              dep_x = "xt", indep_x = "xp", alpha = 0.05, quantile = "t")

sa <- add_bound(sa, arrow = "UD", kind = "direct", lb = -0.3, ub = 0.3)
sa <- add_bound(sa, arrow = "ZY", kind = "direct", lb = -0.3, ub = 0.3)
sa <- add_bound(sa, arrow = "UD", kind = "comparative",
                b = 0.2, J = "xp", I = NULL, name = "bud-comp")
sa <- add_bound(sa, arrow = "UY", kind = "comparative-d",
                b = 0.3, J = "xp", I = NULL, name = "buy-comp")

pl <- b_contours(sa, pir_lower = TRUE,
                 bound1 = "bud-comp", range1 = c(0.05, 0.4),
                 bound2 = "buy-comp", range2 = c(0.1, 0.6),
                 val_interest = 0)
# print(pl)

}
