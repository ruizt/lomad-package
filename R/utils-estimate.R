# Internal utilities for estimate_*() functions
#
# .resid_acov_ar1()   — ACVF of eta - eta^(h) for AR(1) eta
# .resid_vratio_ar1() — its V(2)/V(1), invertible for phi
# .variogram_ar1()    — AR(1) estimation from lag-1/2 variogram


# ACVF of the detrending residual r_t = eta_t - eta_t^(h), for AR(1) eta with
# unit innovation variance. Subtracting the trailing MA makes r an FIR filter
# of eta with weights a = (1 - 1/h, -1/h, ..., -1/h), so
#   gamma_r(l) = sum_{j,k} a_j a_k gamma_eta(l + k - j).
.resid_acov_ar1 <- function(phi, h, lag_max) {
  a <- c(1 - 1 / h, rep(-1 / h, h - 1L))
  g <- function(l) phi^abs(l) / (1 - phi^2)
  vapply(0:lag_max, function(l) {
    idx <- outer(seq_along(a), seq_along(a), function(j, k) l + (k - j))
    sum(outer(a, a) * matrix(g(idx), nrow = length(a)))
  }, numeric(1))
}

# V(2)/V(1) of that residual. Free of sigma^2, and monotone increasing in phi,
# so it can be inverted for phi by root-finding.
.resid_vratio_ar1 <- function(phi, h) {
  g <- .resid_acov_ar1(phi, h, 2L)
  (g[1L] - g[3L]) / (g[1L] - g[2L])
}


# AR(1) estimation from the lag-1 and lag-2 variogram.
#
# `h` is the width of the moving average that was subtracted to form `resid`.
# Leaving it NULL treats `resid` as the noise itself and inverts the plain
# AR(1) relation V(2)/V(1) = 1 + phi. Supplying it corrects for the fact that
# subtracting a trailing MA removes part of the noise along with the trend:
# the leftover eta - eta^(h) is a high-pass-filtered AR(1), not an AR(1), and
# the naive inversion reads its phi too high -- by about +0.08 at h = 5, worth
# 20-30% in the long-run variance L that V depends on.
.variogram_ar1 <- function(resid, h = NULL) {
  n  <- length(resid)
  d1 <- diff(resid)
  d2 <- resid[3:n] - resid[1:(n - 2L)]

  V1 <- mean(d1^2) / 2
  V2 <- mean(d2^2) / 2
  ratio <- V2 / V1

  use_corr <- !is.null(h) && is.finite(h) && h > 1L
  if (use_corr) {
    lo <- .resid_vratio_ar1(0.01, h)
    hi <- .resid_vratio_ar1(0.99, h)
    if (ratio <= lo) {
      phi_raw <- 0.01
    } else if (ratio >= hi) {
      phi_raw <- 0.99
    } else {
      phi_raw <- stats::uniroot(
        function(p) .resid_vratio_ar1(p, h) - ratio,
        interval = c(0.01, 0.99), tol = 1e-8)$root
    }
    boundary <- ratio <= lo || ratio >= hi
  } else {
    phi_raw  <- ratio - 1
    boundary <- phi_raw > 0.99 || phi_raw < 0
  }
  phi <- min(max(phi_raw, 0.01), 0.99)

  # A raw estimate at or beyond the clamp boundaries signals that the AR(1)
  # model does not describe the residual autocovariance (e.g. oscillatory or
  # strongly persistent noise). Silently clamping would return arbitrarily
  # wrong noise variances, so surface it.
  if (boundary) {
    warning(sprintf(paste0(
      "Variogram AR(1) estimate hit its boundary (raw phi = %.3f, clamped ",
      "to %.2f); the AR(1) noise model is likely misspecified for this ",
      "series. Estimate the noise autocovariance yourself and pass it via ",
      "noise_override; see vignette(\"lomad\")."), phi_raw, phi))
  }

  # sigma^2 from V1: on the corrected path V1 is the filtered residual's
  # variogram, so divide out the filter's effect rather than reading V2.
  sigma2 <- if (use_corr) {
    g <- .resid_acov_ar1(phi, h, 1L)
    V1 / (g[1L] - g[2L])
  } else {
    V2
  }

  list(
    ar     = phi,
    ma     = numeric(0),
    sigma2 = sigma2,
    order  = c(1L, 0L, 0L)
  )
}
