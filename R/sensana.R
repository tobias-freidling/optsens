## method that generates the sensana object
## don't supply iv or lm model, just the data
## helper function to compute IV and regression estimates
## plus confidence intervals
## return sensana object

## required y, d, z vectors; x data.frame
sensana <- function(y, d, indep_x, dep_x, x = NULL, z = NULL,
                    intercept = TRUE, alpha = 0.05) {
  cl <- match.call()

  ## LM estimates
  lm_df <- data.frame(y = y, d = d)
  if (!is.null(x)) lm_df <- cbind(lm_df, x)
  if (!is.null(z)) lm_df <- cbind(lm_df, z)
  formula <- if(intercept) "y ~ ." else "y ~ -1 + ."

  model_ols <- lm(formula, data = lm_df)
  beta_ols <- as.vector(model_ols$coefficients["d"])
  confint_ols <- as.vector(confint(model_ols, "d", level = 1 - alpha))

  ## IV estimates
  if (!is.null(z)) {
    model_iv <- ivmodel::ivmodel(Y = y, D = d, X = x, Z = as.matrix(z),
                                 intercept = intercept, alpha = alpha, k = 1)
    kclass <- ivmodel::KClass(model_iv, k = 1, alpha = alpha)
    beta_iv <- as.vector(kclass$point.est)
    confint_iv <- as.vector(kclass$ci)
  } else {
    model_iv <- NULL
    beta_iv <- NULL
    confint_iv <- NULL
  }

  ## Centring and disentangling covariates
  yc <- as.vector(scale(y, scale = FALSE))
  dc <- as.vector(scale(d, scale = FALSE))
  zc <- if(is.null(z)) NULL else as.vector(scale(z, scale = FALSE))
  if (is.null(x) || is.null(dep_x)) {
    xtc <- NULL
  } else {
    xtc <- as.matrix(sapply(x[, dep_x],
                            function(a) scale(a, scale = FALSE)))
  }
  if (is.null(x) || is.null(indep_x)) {
    xpc <- NULL
  } else {
    xpc <- as.matrix(sapply(x[, indep_x],
                            function(a) scale(a, scale = FALSE)))
  }


  ## Initialization of bounds data frame
  bounds <- data.frame(matrix(nrow = 0, ncol = 7))
  colnames(bounds) <- c("arrow", "kind", "lb", "ub", "b", "I", "J")

  sa <- list(call = cl,
             y = yc,
             d = dc,
             xt = xtc,
             xp = xpc,
             z = zc,
             beta_ols = beta_ols,
             beta_iv = beta_iv,
             confint_ols = confint_ols,
             confint_iv = confint_iv,
             model_ols = model_ols,
             model_iv = model_iv,
             bounds = bounds,
             rbc = 0 ## running bound counter
             )
  class(sa) <- "sensana"
  sa
}
