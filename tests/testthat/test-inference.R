
test_that("pir(...) works", {
  set.seed(2024)
  u <- rnorm(20)
  xt <- -0.25*u + rnorm(20)
  xp <- 0.5 * xt + rnorm(20)
  z <- 0.25 * xt - xp + u + rnorm(20)
  d <- 0.5 * xt + 0.5 * xp + u + 2*z + rnorm(20)
  y <- d + 2*xp - xt + z + u + rnorm(20)
  
  xtc <- xt - mean(xt)
  xpc <- xp - mean(xp)
  zc <- z - mean(z)
  dc <- d - mean(d)
  yc <- y - mean(y)
  
  sa <- sensana(y = y, d = d, x = data.frame(xt = xt, xp = xp), z = z,
                dep_x = "xt", indep_x = "xp", alpha = 0.05,
                quantile = "t")
  
  ## Direct bounds
  sa <- add_bound(sa, arrow = "UD", kind = "direct", lb = -0.3, ub = 0.3)
  sa <- add_bound(sa, arrow = "UY", kind = "direct", lb = -0.5, ub = 0.5)
  sa <- add_bound(sa, arrow = "ZU", kind = "direct", lb = -0.3, ub = 0.3)
  sa <- add_bound(sa, arrow = "ZY", kind = "direct", lb = -0.6, ub = 0.6)
  
  ## Indirect bounds
  sa <- add_bound(sa, arrow = "UD", kind = "comparative",
                  b = 0.7 * (1 - r2(dc,xpc,cbind(xtc,zc)) ) / r2(dc,xpc,cbind(xtc,zc)),
                  J = "xp")
  sa <- add_bound(sa, arrow = "UY", kind = "comparative",
                  b = 0.6 * (1 - r2(yc,xpc,cbind(xtc,zc)) ) / r2(yc,xpc,cbind(xtc,zc)), J = "xp")
  sa <- add_bound(sa, arrow = "UY", kind = "comparative-d", J = "xp",
                  b = 0.8 / r2(yc, xpc, cbind(dc,xtc,zc)))
  sa <- add_bound(sa, arrow = "ZU", kind = "comparative", J = "xp",
                  b = 0.5 / r2(zc, xpc, xtc))
  
  pir_sa <- pir(sa)
  
  expect_type(pir_sa, "double")
  expect_true(is.vector(pir_sa))
  expect_equal(length(pir_sa), 2)
})



test_that("Computing collapsing PIR correctly", {
  set.seed(2024)
  u <- rnorm(20)
  x <- rnorm(20)
  d <- x + u
  y <- d + 2*x + u
  
  xc <- x - mean(x)
  dc <- d - mean(d)
  yc <- y - mean(y)
  
  sa <- sensana(y = y, d = d, x = data.frame(x = x), z = NULL,
                dep_x = NULL, indep_x = "x", alpha = 0.05,
                quantile = "t")
  
  sa <- add_bound(sa, arrow = "UD", kind = "direct", lb = -0.3, ub = 0.3)
  sa <- add_bound(sa, arrow = "UY", kind = "direct", lb = -0.5, ub = 0.5)
  sa <- add_bound(sa, arrow = "UD", kind = "comparative",
                  b = 0.7 * (1-r2(dc,xc)) / r2(dc,xc), J = "x")
  sa <- add_bound(sa, arrow = "UY", kind = "comparative",
                  b = 0.6 * (1-r2(yc,xc)) / r2(yc,xc), J = "x")
  sa <- add_bound(sa, arrow = "UY", kind = "comparative-d",
                  b = 0.8 / r2(yc, xc, dc), J = "x")
  
  pir_sa <- pir(sa)
  
  expect_equal(pir_sa, c(2,2))
})



test_that("sensint(...) works", {
  set.seed(2024)
  u <- rnorm(20)
  xt <- -0.25*u + rnorm(20)
  xp <- 0.5 * xt + rnorm(20)
  z <- 0.25 * xt - xp + u + rnorm(20)
  d <- 0.5 * xt + 0.5 * xp + u + 2*z + rnorm(20)
  y <- d + 2*xp - xt + z + u + rnorm(20)
  
  xtc <- xt - mean(xt)
  xpc <- xp - mean(xp)
  zc <- z - mean(z)
  dc <- d - mean(d)
  yc <- y - mean(y)
  
  sa <- sensana(y = y, d = d, x = data.frame(xt = xt, xp = xp), z = z,
                dep_x = "xt", indep_x = "xp", alpha = 0.05,
                quantile = "t")
  
  ## Direct bounds
  sa <- add_bound(sa, arrow = "UD", kind = "direct", lb = -0.3, ub = 0.3)
  sa <- add_bound(sa, arrow = "UY", kind = "direct", lb = -0.5, ub = 0.5)
  sa <- add_bound(sa, arrow = "ZU", kind = "direct", lb = -0.3, ub = 0.3)
  sa <- add_bound(sa, arrow = "ZY", kind = "direct", lb = -0.6, ub = 0.6)
  
  ## Indirect bounds
  sa <- add_bound(sa, arrow = "UD", kind = "comparative",
                  b = 0.7 * (1 - r2(dc,xpc,cbind(xtc,zc)) ) / r2(dc,xpc,cbind(xtc,zc)),
                  J = "xp")
  sa <- add_bound(sa, arrow = "UY", kind = "comparative",
                  b = 0.6 * (1 - r2(yc,xpc,cbind(xtc,zc)) ) / r2(yc,xpc,cbind(xtc,zc)), J = "xp")
  sa <- add_bound(sa, arrow = "UY", kind = "comparative-d", J = "xp",
                  b = 0.8 / r2(yc, xpc, cbind(dc,xtc,zc)))
  sa <- add_bound(sa, arrow = "ZU", kind = "comparative", J = "xp",
                  b = 0.5 / r2(zc, xpc, xtc))
  
  sensint_sa <- sensint(sa, alpha = 0.2,
                        boot_procedure = c("norm", "basic", "perc", "bca"),
                        boot_samples = 50)
  
  expect_s3_class(sensint_sa, "sensint")
  expect_s3_class(sensint_sa$boot_obj, "boot")
  expect_equal(sensint_sa$alpha, 0.2)
  expect_true(is.data.frame(sensint_sa$sensint))
  expect_type(sensint_sa$sensint$sl, "double")
  expect_type(sensint_sa$sensint$su, "double")
})

