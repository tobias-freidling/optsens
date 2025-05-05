


test_that("sensana object is correct", {
  x <- rnorm(20)
  z <- 0.5 * x + rnorm(20)
  d <- z - 0.5 * x + rnorm(20)
  y <- d + x + rnorm(20)
  
  sa <- sensana(y = y, d = d, x = data.frame(x = x), z = z,
                dep_x = NULL, indep_x = "x", alpha = 0.05,
                quantile = "t")
  
  expect_s3_class(sa, "sensana")
  expect_true(is.data.frame(sa$bounds))
  expect_true(is.vector(sa$y))
  expect_true(is.vector(sa$d))
  expect_true(is.vector(sa$z))
  expect_true(is.matrix(sa$xp))
  expect_true(is.null(sa$xt))
})




test_that("sensana inference is correct", {
  x <- rnorm(20)
  z <- 0.5 * x + rnorm(20)
  d <- z - 0.5 * x + rnorm(20)
  y <- d + x + rnorm(20)
  
  ## sensana
  sa <- sensana(y = y, d = d, x = data.frame(x = x), z = z,
                dep_x = NULL, indep_x = "x", alpha = 0.05,
                quantile = "t")
  
  ## OLS inference
  model_ols <- stats::lm("y ~ x + z + d",
                         data = data.frame(y = y, x = x, d = d, z = z))
  beta_lm <- as.vector(model_ols$coefficients["d"])
  confint_lm <- as.vector(stats::confint(model_ols, "d", level = 0.95))
  
  
  ## TSLS inference
  model_iv <- ivmodel::ivmodel(Y = y, D = d, X = x, Z = as.matrix(z),
                               intercept = TRUE, alpha = 0.05, k = 1)
  kclass <- ivmodel::KClass(model_iv, k = 1, alpha = 0.05)
  beta_iv <- as.vector(kclass$point.est)
  confint_iv <- as.vector(kclass$ci)
  
  expect_equal(sa$beta_ols, beta_lm)
  expect_equal(sa$confint_ols, confint_lm)
  expect_equal(sa$beta_tsls, beta_iv)
  expect_equal(sa$confint_tsls, confint_iv)
})
