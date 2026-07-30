# Internal utilities for estimate_*() functions
#
# .variogram_ar1()    — AR(1) estimation from lag-1/2 variogram


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
      "series. Estimate the noise autocovariance yourself and pass it via ",
      "noise_override; see vignette(\"lomad\")."), phi_raw, phi))
  }

  list(
    ar     = phi,
    ma     = numeric(0),
    sigma2 = V2,
    order  = c(1L, 0L, 0L)
  )
}

