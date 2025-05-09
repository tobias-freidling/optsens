
<!-- README.md is generated from README.Rmd. Please edit that file -->

# optsens

<!-- badges: start -->

<!-- badges: end -->

This package provides tools to conduct sensitivity analysis for the
identification assumptions of linear regression and instrumental
variable models. Users provide the data and specify direct or
comparative bounds on the influence of an unmeasured confounder on
observed variables. The package estimates the resulting partially
identified region (PIR) and provides uncertainty quantification.
Moreover, the results can be visualized via 2 types of contour plots (R-
and b-contours).

This is accomplished by casting sensitivity analysis as a constraint
stochastic optimization problem; optsens uses a tailored grid-search
algorithm to solve the resulting optimization problem.

optsens is based on the article *Optimization-based Sensitivity Analysis
for Unmeasured Confounding using Partial Correlations*, available on
[arXiv](https://arxiv.org/abs/2301.00040). The code to create the plots
and simulation studies in the paper can be found in this [Github
repository](https://github.com/tobias-freidling/optsens-replication).

## Installation

You can install the development version of optsens from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tobias-freidling/optsens")
```
