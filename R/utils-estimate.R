# Internal utilities for estimate_*() functions
#
# .variogram_ar1()    — AR(1) estimation from lag-1/2 variogram
# .yule_walker()      — solve Yule-Walker equations


# AR(1) estimation from the lag-1 and lag-2 variogram.
.variogram_ar1 <- function(resid) {
  n  <- length(resid)
  d1 <- diff(resid)
  d2 <- resid[3:n] - resid[1:(n - 2L)]

  V1 <- mean(d1^2) / 2
  V2 <- mean(d2^2) / 2

  phi_raw <- V2 / V1 - 1
  phi     <- min(max(phi_raw, 0.01), 0.99)

  # A raw estimate at or beyond the clamp boundaries signals that the AR(1)
  # model does not describe the residual autocovariance (e.g. oscillatory or
  # strongly persistent noise). Silently clamping would return arbitrarily
  # wrong noise variances, so surface it.
  if (phi_raw > 0.99 || phi_raw < 0) {
    warning(sprintf(paste0(
      "Variogram AR(1) estimate hit its boundary (raw phi = %.3f, clamped ",
      "to %.2f); the AR(1) noise model is likely misspecified for this ",
      "series. Consider noise_method = \"arma\" or supplying ",
      "noise_override."), phi_raw, phi))
  }

  list(
    ar     = phi,
    ma     = numeric(0),
    sigma2 = V2,
    order  = c(1L, 0L, 0L)
  )
}


# Solve Yule-Walker equations for AR(p).
.yule_walker <- function(gamma0, gamma_hat, p) {
  acvf <- c(gamma0, gamma_hat[seq_len(p)])
  Gamma_mat <- stats::toeplitz(acvf[seq_len(p)])
  gamma_vec <- gamma_hat[seq_len(p)]

  phi <- tryCatch(
    solve(Gamma_mat, gamma_vec),
    error = function(e) rep(0, p)
  )

  sigma2 <- gamma0 - sum(phi * gamma_vec)
  sigma2 <- max(sigma2, 1e-10)

  list(ar = as.numeric(phi), sigma2 = sigma2)
}
