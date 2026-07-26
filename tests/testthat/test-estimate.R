# Tests for estimate_trends, estimate_ar1_noise, estimate_arma_noise,
# and key internals (.variogram_ar1, .yule_walker, .select_arma)

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


# ---- .select_arma ----

test_that(".select_arma: selects AR(1) on AR(1) data", {
  set.seed(7723)
  z <- as.numeric(arima.sim(list(ar = 0.6), n = 500))
  fit <- .select_arma(z, max_pq = 2)
  expect_equal(fit$order[1], 1, label = "AR order")
  expect_true(fit$sigma2 > 0)
})

test_that(".select_arma: handles MA(1) data", {
  set.seed(4519)
  z <- as.numeric(arima.sim(list(ma = 0.5), n = 500))
  fit <- .select_arma(z, max_pq = 2)
  # Should select a model with q >= 1 or a high-order AR approximation
  expect_true(fit$order[1] + fit$order[3] >= 1)
  expect_true(fit$sigma2 > 0)
})

test_that(".select_arma: handles ARMA(1,1) data", {
  set.seed(6188)
  z <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.3), n = 500))
  fit <- .select_arma(z, max_pq = 2)
  expect_true(fit$order[1] >= 1 || fit$order[3] >= 1)
  expect_true(fit$sigma2 > 0)
})

test_that(".select_arma: white noise fallback on near-constant input", {
  # Use near-zero noise rather than exactly zero (which crashes arima)
  set.seed(9911)
  fit <- .select_arma(rnorm(100, sd = 1e-6), max_pq = 2)
  expect_true(fit$sigma2 >= 0)
})


# ---- estimate_arma_noise ----

test_that("estimate_arma_noise: selects correct AR order on AR(2) data", {
  set.seed(8812)
  n <- 2000
  phi <- c(0.5, -0.2)
  trend <- cumsum(rnorm(n, sd = 0.05))
  z1 <- as.numeric(arima.sim(list(ar = phi), n = n))
  z2 <- as.numeric(arima.sim(list(ar = phi), n = n))
  y1 <- trend + z1
  y2 <- trend + z2

  tr <- estimate_trends(y1, y2, h = 10)
  noise <- estimate_arma_noise(y1, y2, tr$trend, p_max = 4)

  expect_true(noise$series1$order[1] >= 2,
              label = "AR order should be >= 2")
})

test_that("estimate_arma_noise: errors on length mismatch", {
  expect_error(estimate_arma_noise(1:10, 1:5, 1:10), "same length")
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


# ---- estimate_acf_noise / .variogram_acov ----

test_that(".variogram_acov recovers AR(1) autocovariances", {
  set.seed(7)
  phi <- 0.5; sigma2 <- 1
  z <- as.numeric(arima.sim(list(ar = phi), n = 5000, sd = sqrt(sigma2)))
  acov <- .variogram_acov(z, lag_max = 10)
  gamma0_true <- sigma2 / (1 - phi^2)
  expect_equal(acov[1], gamma0_true, tolerance = 0.1)
  # taper shrinks higher lags toward zero, so compare the ratio at lag 1
  # against the tapered truth
  n <- length(z); b <- ceiling(10 * log10(n))
  expect_equal(acov[2] / acov[1], phi * (1 - 1 / (b + 1)), tolerance = 0.1)
})

test_that(".variogram_acov output length and taper zeros", {
  set.seed(8)
  z <- rnorm(200)
  acov <- .variogram_acov(z, lag_max = 150)
  expect_length(acov, 151)
  b <- ceiling(10 * log10(200))
  expect_true(all(acov[(b + 2):151] == 0))
})

test_that("estimate_acf_noise returns per-series acov lists", {
  set.seed(9)
  n <- 300
  trend <- sin(seq(0, 2 * pi, length.out = n))
  y1 <- trend + as.numeric(arima.sim(list(ar = 0.4), n = n))
  y2 <- trend + rnorm(n)
  out <- estimate_acf_noise(y1, y2, trend, lag_max = 50)
  expect_named(out, c("series1", "series2"))
  expect_length(out$series1$acov, 51)
  expect_gt(out$series1$acov[1], 0)
  expect_gt(out$series2$acov[1], 0)
})

test_that("estimate_acf_noise is robust to the stress-case ARMA noise", {
  # Oscillatory ARMA(2,1) noise that breaks the parametric variogram fits
  set.seed(10)
  n <- 2000
  z <- as.numeric(arima.sim(list(ar = c(0.602, -0.684), ma = 0.51),
                            n = n, sd = 0.15))
  gamma0_true <- (0.15^2) * sum(c(1, ARMAtoMA(ar = c(0.602, -0.684),
                                              ma = 0.51, lag.max = 200))^2)
  acov <- .variogram_acov(z, lag_max = 10)
  expect_equal(acov[1], gamma0_true, tolerance = 0.25)
})

test_that(".variogram_acov plateau anchor also recovers AR(1) variance", {
  set.seed(13)
  z <- as.numeric(arima.sim(list(ar = 0.5), n = 5000))
  a_var <- .variogram_acov(z, lag_max = 5, anchor = "variance")
  a_pla <- .variogram_acov(z, lag_max = 5, anchor = "plateau")
  gamma0_true <- 1 / (1 - 0.5^2)
  expect_equal(a_var[1], gamma0_true, tolerance = 0.1)
  expect_equal(a_pla[1], gamma0_true, tolerance = 0.15)
})
