# Tests for estimate_trends, estimate_ar1_noise,
# and key internals (.variogram_ar1)

# ---- estimate_trends ----

test_that("estimate_trends: returns correct structure", {
  set.seed(3821)
  n <- 200
  trend <- cumsum(rnorm(n, sd = 0.1))
  y1 <- trend + rnorm(n, sd = 0.5)
  y2 <- trend + rnorm(n, sd = 0.5)

  tr <- estimate_trends(y1, y2, h = 10)
  expect_true(all(c("ma1", "ma2") %in% names(tr)))
  expect_length(tr$ma1, n)
  expect_true(any(is.na(tr$ma1)))
})


# ---- estimate_ar1_noise ----

test_that("estimate_ar1_noise: recovers phi approximately", {
  set.seed(5591)
  n <- 2000
  true_phi <- 0.6
  sigma2 <- 1
  trend <- cumsum(rnorm(n, sd = 0.05))
  z1 <- as.numeric(arima.sim(list(ar = true_phi), n = n, sd = sqrt(sigma2)))
  z2 <- as.numeric(arima.sim(list(ar = true_phi), n = n, sd = sqrt(sigma2)))
  y1 <- trend + z1
  y2 <- trend + z2

  tr <- estimate_trends(y1, y2, h = 10)
  noise <- estimate_ar1_noise(y1, y2)

  expect_equal(noise$series1$ar, true_phi, tolerance = 0.15,
               label = "AR(1) phi recovery")
  expect_true(noise$series1$sigma2 > 0)
  expect_equal(noise$series1$order, c(1, 0, 0))
})


# ---- .variogram_ar1 ----

test_that(".variogram_ar1: known AR(1) within tolerance", {
  set.seed(2781)
  true_phi <- 0.5
  z <- as.numeric(arima.sim(list(ar = true_phi), n = 3000))
  est <- .variogram_ar1(z)
  expect_equal(est$ar, true_phi, tolerance = 0.15)
})



# ---- .variogram_ar1 boundary warning ----

test_that(".variogram_ar1 warns when the AR(1) fit hits a clamp boundary", {
  set.seed(42)
  # Strongly oscillatory residuals: lag-1 differences are large relative to
  # lag-2 differences, driving the raw phi estimate negative.
  resid_osc <- rep(c(1, -1), 60) + rnorm(120, sd = 0.05)
  expect_warning(.variogram_ar1(resid_osc), "boundary")

  # Well-behaved AR(1) residuals: no warning.
  set.seed(43)
  resid_ok <- as.numeric(arima.sim(list(ar = 0.5), n = 400))
  expect_no_warning(.variogram_ar1(resid_ok))
})
