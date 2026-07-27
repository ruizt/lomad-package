# Tests for estimate_trends, estimate_ar1_noise,
# and key internals (.variogram_ar1, .yule_walker)

# ---- estimate_trends ----

test_that("estimate_trends: returns correct structure", {
  set.seed(3821)
  n <- 200
  trend <- cumsum(rnorm(n, sd = 0.1))
  y1 <- trend + rnorm(n, sd = 0.5)
  y2 <- trend + rnorm(n, sd = 0.5)

  tr <- estimate_trends(y1, y2, h = 10)
  expect_true(all(c("trend", "ma1", "ma2") %in% names(tr)))
  expect_length(tr$trend, n)
  expect_true(any(is.na(tr$trend)))
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
  noise <- estimate_ar1_noise(y1, y2, tr$trend)

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


# ---- .yule_walker ----

test_that(".yule_walker: exact solution on known Toeplitz system", {
  # AR(2) with known coefficients
  phi <- c(0.5, -0.3)
  sigma2 <- 1
  # gamma(0) from AR(2): gamma(0) = sigma2 / (1 - phi1*rho1 - phi2*rho2)
  # Use ARMAacf to get the exact ACF
  acf_vals <- as.numeric(ARMAacf(ar = phi, lag.max = 5))
  gamma0 <- sigma2 / (1 - sum(phi * acf_vals[2:3]))
  gamma_hat <- gamma0 * acf_vals[2:6]

  yw <- .yule_walker(gamma0, gamma_hat, p = 2)
  expect_equal(yw$ar, phi, tolerance = 1e-6)
  expect_equal(yw$sigma2, sigma2, tolerance = 1e-4)
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
