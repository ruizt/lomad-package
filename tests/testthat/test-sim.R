# Tests for sim_trends(), sim_noise(), sim_noise_pair(), and key internals

# ---- sim_trends ----

test_that("sim_trends: returns correct structure and length", {
  tr <- sim_trends(n = 200, d = 1, seed = 4831)
  expect_type(tr, "list")
  expect_named(tr, c("x1", "x2", "x_mean", "w", "a", "b"))
  expect_length(tr$x1, 200)
  expect_length(tr$w, 200)
  expect_equal(tr$a, rep(0, 200))
  expect_equal(tr$b, rep(1, 200))
})

test_that("sim_trends: separation is linear in d, not normalised to it", {
  # d scales the distinct component; it is deliberately not solved for so that
  # ||x1 - x2|| == d, which would amplify structures that decouple briefly.
  for (m in c("dist", "rs", "rm", "fr")) {
    args <- list(200, method = m, seed = 7713)
    args <- c(args, if (m == "fr") list(rate = 0.02) else
                    if (m == "dist") NULL else list(bw = 30, coupling = 0.8))
    set.seed(2); a <- do.call(sim_trends, c(args, list(d = 1)))
    set.seed(2); b <- do.call(sim_trends, c(args, list(d = 3)))
    expect_equal(b$x1 - b$x2, 3 * (a$x1 - a$x2), tolerance = 1e-10,
                 label = sprintf("linearity, method = %s", m))
    # and the coupling profile is untouched by d
    expect_equal(a$w, b$w, label = sprintf("w invariance, method = %s", m))
  }
})

test_that("sim_trends: d = 0 collapses both series onto the shared mean", {
  tr <- sim_trends(200, d = 0, method = "rs", bw = 30, seed = 7713)
  expect_equal(tr$x1, tr$x_mean)
  expect_equal(tr$x2, tr$x_mean)
})

test_that("sim_trends: method='dist' gives w = 0", {
  tr <- sim_trends(100, d = 1, method = "dist", seed = 5122)
  expect_true(all(tr$w == 0))
})

test_that("sim_trends: method='smooth' gives w in (0, 1)", {
  tr <- sim_trends(500, d = 1, method = "rs", bw = 50,
                   coupling = 0.8, seed = 5122)
  expect_true(all(tr$w > 0 & tr$w < 1))
})

test_that("sim_trends: method='rate' gives w in [0, 1]", {
  tr <- sim_trends(500, d = 1, method = "fr", rate = 0.01, seed = 5122)
  expect_true(all(tr$w >= 0 & tr$w <= 1))
})

test_that("sim_trends: custom w (vector and function)", {
  w_vec <- rep(c(0, 1), each = 100)
  tr_v <- sim_trends(200, d = 2, w = w_vec, seed = 3481)
  expect_equal(tr_v$w, w_vec)
  # w = 1 couples the pair exactly, w = 0 leaves the full distinct component
  expect_equal(tr_v$x1[101:200], tr_v$x2[101:200])
  expect_true(all(tr_v$x1[1:100] != tr_v$x2[1:100]))

  tr_f <- sim_trends(200, d = 2, w = \(n) rep(0.5, n), seed = 3481)
  expect_equal(tr_f$w, rep(0.5, 200))
})

test_that("sim_trends: errors on bad inputs", {
  expect_error(sim_trends(100, nb = 24), "`nb` must be odd")
  expect_error(sim_trends(100, w = rep(0, 50)), "`w` must have length")
})

test_that("sim_trends: seed reproducibility", {
  tr1 <- sim_trends(100, d = 1, seed = 8192)
  tr2 <- sim_trends(100, d = 1, seed = 8192)
  expect_identical(tr1, tr2)
})


# ---- .generate_coef_pair ----

test_that(".generate_coef_pair: displacement is orthogonal to coef1", {
  for (d in c(0.5, 1, 3)) {
    cp <- .generate_coef_pair(nb = 25, d = d, seed = 2244)
    # coef2 is rescaled, so recover the pre-rescaling displacement direction
    # by comparing against the unit vector along coef1.
    u1  <- cp$coef1 / sqrt(sum(cp$coef1^2))
    off <- cp$coef2 - sum(cp$coef2 * u1) * u1     # component orthogonal to c1
    expect_gt(sqrt(sum(off^2)), 0)
  }
})

test_that(".generate_coef_pair: equal norms, so lambda_1 == lambda_2", {
  for (d in c(0, 0.5, 1, 3)) {
    cp <- .generate_coef_pair(nb = 25, d = d, seed = 2244)
    expect_equal(sqrt(sum(cp$coef2^2)), sqrt(sum(cp$coef1^2)),
                 tolerance = 1e-8,
                 label = sprintf("coef2 norm, d = %g", d))
  }
})

test_that(".generate_coef_pair: d maps to the affine effect size", {
  # Orthogonal displacement then rescaling preserves the angle, so
  #   r = cos(angle) = ||c1|| / sqrt(||c1||^2 + d^2)
  # exactly, and delta = sqrt(1 - r^2) = d / sqrt(||c1||^2 + d^2).
  for (d in c(0.5, 1, 3)) {
    cp <- .generate_coef_pair(nb = 25, d = d, seed = 2244)
    n1 <- sqrt(sum(cp$coef1^2))
    r  <- sum(cp$coef1 * cp$coef2) / (n1 * sqrt(sum(cp$coef2^2)))
    expect_equal(r, n1 / sqrt(n1^2 + d^2), tolerance = 1e-8,
                 label = sprintf("r, d = %g", d))
  }
})

test_that(".generate_coef_pair: d = 0 is exactly the null", {
  cp <- .generate_coef_pair(nb = 25, d = 0, seed = 2244)
  expect_equal(cp$coef2, cp$coef1, tolerance = 1e-10)
})


# ---- .make_w_smooth ----

test_that(".make_w_smooth: output in (0, 1)", {
  set.seed(6012)
  w <- .make_w_smooth(500, bw = 50, coupling = 0.8)
  expect_length(w, 500)
  expect_true(all(w > 0 & w < 1))
})


# ---- .make_w_rate ----

test_that(".make_w_rate: output in [0, 1]", {
  w <- .make_w_rate(500, rate = 0.01)
  expect_length(w, 500)
  expect_true(all(w >= 0 & w <= 1))
})


# ---- sim_noise / sim_noise_pair ----

test_that("sim_noise_pair: returns expected structure", {
  tr <- sim_trends(200, d = 1, seed = 9312)
  sp <- suppressMessages(
    sim_noise_pair(tr, h = 10, ar.coefs = 0.3, lambda_target = 2, seed = 9312)
  )
  expect_named(sp, c("y1", "y2", "x1", "x2", "noise"))
  expect_length(sp$y1, 200)
  expect_true(sp$noise$series1$sigma > 0)
})

test_that("sim_noise_pair: seed reproducibility", {
  tr <- sim_trends(200, d = 1, seed = 1142)
  sp1 <- suppressMessages(
    sim_noise_pair(tr, h = 10, ar.coefs = 0.3, lambda_target = 1, seed = 4401)
  )
  sp2 <- suppressMessages(
    sim_noise_pair(tr, h = 10, ar.coefs = 0.3, lambda_target = 1, seed = 4401)
  )
  expect_identical(sp1$y1, sp2$y1)
})


# ---- rate coupling weight: pulse shape ----

test_that("`bump` changes pulse shape without moving the events", {
  wg <- lomad:::.make_w_rate(600, rate = 0.01, bump = "gamma")
  wn <- lomad:::.make_w_rate(600, rate = 0.01, bump = "gaussian")
  expect_length(wn, 600)
  expect_true(all(wn >= 0 & wn <= 1))
  # same number of dips, in the same neighbourhoods
  dips <- function(w) which(diff(sign(diff(w))) > 0) + 1L
  expect_equal(length(dips(wg)), length(dips(wn)))
  expect_true(max(abs(sort(dips(wg)) - sort(dips(wn)))) < 20L)
  expect_false(isTRUE(all.equal(wg, wn)))
})

test_that("the gaussian pulse leaks less through differencing", {
  # The point of the option: a corner at onset survives differencing, a smooth
  # pulse does not. Compare each weight's variogram at short lags.
  leak <- function(w) {
    v <- vapply(1:10, function(l)
      mean((w[(l + 1):length(w)] - w[1:(length(w) - l)])^2) / 2, numeric(1))
    v[10] / stats::var(w)
  }
  expect_lt(leak(lomad:::.make_w_rate(600, bump = "gaussian")),
            leak(lomad:::.make_w_rate(600, bump = "gamma")))
})

test_that("sim_trends passes `bump` through and defaults to gaussian", {
  a <- sim_trends(600, d = 1, method = "fr", seed = 7)
  b <- sim_trends(600, d = 1, method = "fr", seed = 7, bump = "gaussian")
  expect_equal(a$x1, b$x1)
  g <- sim_trends(600, d = 1, method = "fr", seed = 7, bump = "gamma")
  expect_false(isTRUE(all.equal(a$x1, g$x1)))
  expect_true(all(is.finite(g$x1)) && all(is.finite(g$x2)))
})


# ---- affine layer ----

test_that("sim_trends: affine layer leaves x1 and the base trends untouched", {
  set.seed(1); a <- sim_trends(400, d = 1, method = "rs", bw = 30, seed = 22)
  set.seed(1); b <- sim_trends(400, d = 1, method = "rs", bw = 30, seed = 22,
                               affine_s = 80)
  expect_identical(a$x1, b$x1)
  expect_identical(a$w, b$w)
  expect_false(identical(a$x2, b$x2))
})

test_that("sim_trends: affine drift respects the per-window cap", {
  tr <- sim_trends(2000, d = 0, method = "rs", bw = 50, seed = 31,
                   affine_s = 100, affine_cap = 0.02)
  expect_equal(max(abs(diff(tr$b, lag = 100))), 0.02, tolerance = 1e-8)
  expect_equal(max(abs(diff(tr$a, lag = 100))), 0.02 * sd(tr$x1),
               tolerance = 1e-8)
  expect_gt(diff(range(tr$b)), 0.02)   # accumulates beyond one window
})

test_that("sim_trends: affine arguments are validated", {
  expect_error(sim_trends(200, affine_s = 1), "between 2")
  expect_error(sim_trends(200, affine_s = 500), "between 2")
  expect_error(sim_trends(200, affine_s = 50, affine_bw = 0), "positive")
  expect_error(sim_trends(200, affine_s = 50, affine_cap = -1), "non-negative")
})

# ---- .compute_delta ----

test_that(".compute_delta: least-squares branch is sqrt(1 - r^2)", {
  set.seed(9)
  x1 <- cumsum(rnorm(400))
  x2 <- 0.4 + 1.3 * x1 + 0.3 * cumsum(rnorm(400))
  s  <- 60L
  expected <- vapply(seq_along(x1), function(t) {
    if (t < s) return(NA_real_)
    w <- (t - s + 1L):t
    sqrt(1 - stats::cor(x1[w], x2[w])^2)
  }, numeric(1))
  expect_equal(.compute_delta(x1, x2, s), expected)
})

test_that(".compute_delta: an exact affine map gives zero separation", {
  set.seed(10)
  mu <- cumsum(rnorm(400))
  y  <- -2 + 3 * mu
  expect_equal(max(.compute_delta(mu, y, 60L), na.rm = TRUE), 0,
               tolerance = 1e-12)
  expect_equal(max(.compute_delta(mu, y, 60L, b = rep(3, 400)), na.rm = TRUE), 0,
               tolerance = 1e-12)
})

test_that(".compute_delta: least-squares branch is affine invariant", {
  set.seed(11)
  x1 <- cumsum(rnorm(300)); x2 <- cumsum(rnorm(300))
  expect_equal(.compute_delta(x1, x2, 50L),
               .compute_delta(5 * x1 - 2, -3 * x2 + 9, 50L))
})

test_that(".compute_delta: window convention matches compute_tau_sq", {
  set.seed(12)
  x1 <- cumsum(rnorm(200)); x2 <- cumsum(rnorm(200))
  d  <- .compute_delta(x1, x2, 40L)
  expect_true(all(is.na(d[1:39])))
  expect_true(all(is.finite(d[40:200])))
  expect_length(d, 200)
})

test_that("sim_trends: supplied coefficients recover a fixed pre-existing map", {
  # The validation study builds its pair as trend2 <- 2 * tr$x1 outside the
  # package. Supplying the coefficients must reproduce that bit for bit, which
  # requires drawing nothing.
  n <- 400
  base <- sim_trends(n, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381)
  got  <- sim_trends(n, d = 0, nb = 25, sd0 = 50, p = 1.5, seed = 5381,
                     affine_a = rep(0, n), affine_b = rep(2, n))
  expect_identical(got$x1, base$x1)
  expect_equal(got$x2, 2 * base$x1)
  expect_equal(got$b, rep(2, n))
})

test_that("sim_trends: supplied coefficients consume no RNG", {
  n <- 300
  set.seed(77); a <- sim_trends(n, d = 1, method = "rs", bw = 30, seed = 5)
  set.seed(77); b <- sim_trends(n, d = 1, method = "rs", bw = 30, seed = 5,
                                affine_b = rep(3, n))
  expect_identical(a$x1, b$x1)
  expect_equal(b$x2, 3 * a$x2)
  # the RNG stream is left in the same place either way
  set.seed(77); invisible(sim_trends(n, d = 1, method = "rs", bw = 30, seed = 5))
  r1 <- runif(1)
  set.seed(77); invisible(sim_trends(n, d = 1, method = "rs", bw = 30, seed = 5,
                                     affine_b = rep(3, n)))
  expect_identical(runif(1), r1)
})

test_that("sim_trends: supplied coefficients override affine_s", {
  n <- 300
  tr <- sim_trends(n, d = 0, seed = 8, affine_s = 60, affine_b = rep(2, n))
  expect_equal(tr$b, rep(2, n))
  expect_equal(tr$a, rep(0, n))
})

test_that("sim_trends: affine coefficients accept functions and validate length", {
  n <- 200
  tr <- sim_trends(n, d = 0, seed = 8, affine_b = function(k) rep(1.5, k))
  expect_equal(tr$b, rep(1.5, n))
  expect_error(sim_trends(n, d = 0, affine_b = rep(2, 10)), "length")
  expect_error(sim_trends(n, d = 0, affine_a = "x"), "numeric")
})

test_that("sim_trends: affine_cap = 0 is a drift-free layer that draws nothing", {
  n <- 300
  set.seed(5); off  <- sim_trends(n, d = 1, method = "rs", bw = 30, seed = 9)
  set.seed(5); zero <- sim_trends(n, d = 1, method = "rs", bw = 30, seed = 9,
                                  affine_s = 60, affine_cap = 0)
  expect_identical(off$x1, zero$x1)
  expect_identical(off$x2, zero$x2)
  expect_equal(zero$a, rep(0, n))
  expect_equal(zero$b, rep(1, n))

  # and the RNG stream is left where the layer-off call leaves it
  set.seed(5); invisible(sim_trends(n, d = 1, method = "rs", bw = 30, seed = 9))
  r1 <- runif(1)
  set.seed(5); invisible(sim_trends(n, d = 1, method = "rs", bw = 30, seed = 9,
                                    affine_s = 60, affine_cap = 0))
  expect_identical(runif(1), r1)
})
