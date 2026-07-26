# Internal utilities for estimate_*() functions
#
# .variogram_ar1()    — AR(1) estimation from lag-1/2 variogram
# .variogram_ar()     — general AR(p) from multi-lag variogram + Yule-Walker
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


# General AR(p) estimation from multi-lag variogram + Yule-Walker.
.variogram_ar <- function(resid, p = NULL, p_max = 5L, L = 20L) {
  n <- length(resid)
  if (L >= n - 1L) L <- n - 2L
  if (is.null(p)) {
    max_order <- p_max
  } else {
    max_order <- p
  }
  if (L < max_order + 1L) L <- max_order + 1L

  V <- numeric(L)
  for (l in seq_len(L)) {
    dl <- resid[(l + 1L):n] - resid[1:(n - l)]
    V[l] <- mean(dl^2) / 2
  }

  plateau_start <- max(1L, ceiling(L / 2))
  gamma0 <- mean(V[plateau_start:L])
  gamma_hat <- gamma0 - V

  if (!is.null(p)) {
    yw <- .yule_walker(gamma0, gamma_hat, p)
    return(list(
      ar     = yw$ar,
      ma     = numeric(0),
      sigma2 = yw$sigma2,
      order  = c(p, 0L, 0L)
    ))
  }

  best_bic <- Inf
  best_p   <- 1L
  best_yw  <- NULL

  for (pp in seq_len(p_max)) {
    yw  <- .yule_walker(gamma0, gamma_hat, pp)
    if (yw$sigma2 <= 0) next
    bic <- n * log(yw$sigma2) + pp * log(n)
    if (bic < best_bic) {
      best_bic <- bic
      best_p   <- pp
      best_yw  <- yw
    }
  }

  if (is.null(best_yw)) {
    best_yw <- .yule_walker(gamma0, gamma_hat, 1L)
    best_p  <- 1L
  }

  list(
    ar     = best_yw$ar,
    ma     = numeric(0),
    sigma2 = best_yw$sigma2,
    order  = c(best_p, 0L, 0L)
  )
}


# Nonparametric ACVF from the multi-lag variogram, Bartlett-tapered.
#
# Same difference-based recovery as .variogram_ar — gamma_hat(l) = gamma_hat(0)
# - V(l) with gamma_hat(0) anchored on the variogram plateau — but keeps the
# nonparametric sequence instead of projecting it onto an AR(p) model. A
# Bartlett taper (weights 1 - l/(b + 1), zero beyond b) stabilises the
# long-run sums L, Q, Q12 computed downstream (Newey-West style); without it,
# raw autocovariance sums can go negative.
#
# resid     : residual series (observed minus shared trend estimate)
# lag_max   : length of the returned ACVF (lags 0..lag_max; zero beyond taper)
# bandwidth : Bartlett taper bandwidth b. NULL uses ceiling(10 * log10(n)).
# anchor    : gamma(0) anchor — "variance" (residual sample variance; tight,
#             test stays calibrated, conservative under trend separation) or
#             "plateau" (variogram plateau average; trend-robust but noisy —
#             anchor error propagates coherently to all lags and can make
#             the one-sided test anti-conservative).
.variogram_acov <- function(resid, lag_max, bandwidth = NULL,
                            anchor = c("variance", "plateau")) {
  anchor <- match.arg(anchor)
  n <- length(resid)
  if (is.null(bandwidth)) bandwidth <- ceiling(10 * log10(n))
  b <- max(1L, min(as.integer(bandwidth), n - 22L))

  # Variogram out to the anchor band beyond the taper bandwidth
  L <- min(b + 20L, n - 2L)
  V <- numeric(L)
  for (l in seq_len(L)) {
    dl <- resid[(l + 1L):n] - resid[1:(n - l)]
    V[l] <- mean(dl^2) / 2
  }

  gamma0 <- if (anchor == "variance") {
    mean((resid - mean(resid))^2)
  } else {
    mean(V[max(b + 1L, ceiling(L / 2)):L])
  }
  if (gamma0 <= 0)
    stop("gamma(0) anchor estimate is non-positive; cannot estimate the ",
         "noise ACVF.")

  # If the variogram is still climbing across the far lags, noise
  # correlation may extend beyond the taper bandwidth.
  plateau <- V[max(b + 1L, ceiling(L / 2)):L]
  half <- length(plateau) %/% 2
  if (half >= 3L) {
    rise <- mean(plateau[(half + 1L):length(plateau)]) / mean(plateau[1:half])
    if (is.finite(rise) && rise > 1.25)
      warning(sprintf(paste0(
        "Variogram is still rising at lag %d (late/early ratio %.2f); ",
        "noise correlation may extend beyond the taper bandwidth. Consider ",
        "a larger `bandwidth` or a different noise method."), L, rise))
  }

  gamma <- gamma0 - V                              # lags 1..L
  wts   <- pmax(0, 1 - seq_len(L) / (b + 1L))     # Bartlett taper
  acov  <- c(gamma0, gamma * wts)                 # lags 0..L

  out <- numeric(lag_max + 1L)
  k   <- min(lag_max + 1L, length(acov))
  out[seq_len(k)] <- acov[seq_len(k)]
  out
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
