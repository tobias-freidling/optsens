
test_that("r works", {
  x <- 1:10
  y <- 2*x + 5
  z <- -3*x + 5*y
  xc <- x - mean(x)
  yc <- y - mean(y)
  zc <- z - mean(z)
  
  expect_equal(r(xc,yc), 1)
  expect_equal(r(xc,-yc), -1)
  expect_equal(r(xc,yc,as.matrix(zc)),1)
})



test_that("r2 works", {
  x <- 1:10
  y <- 2*x + 5
  z <- -3*x + 5*y
  
  xc <- x - mean(x)
  yc <- y - mean(y)
  zc <- z - mean(z)
  
  expect_equal(r2(xc, as.matrix(yc)), 1)
  expect_equal(r2(xc, as.matrix(-yc)), 1)
  expect_equal(r2(xc, as.matrix(yc), as.matrix(zc)), 1)
  expect_equal(r2(xc, as.matrix(yc), as.matrix(zc)),
               (r(xc, yc, as.matrix(zc)))^2)
})



test_that("f and finv work", {
  x <- stats::rnorm(10)
  y <- stats::rnorm(10)
  xc <- x - mean(x)
  yc <- y - mean(y)
  
  r_yx <- r(yc,xc)
  expect_equal(finv(f(r_yx)), r_yx)
})