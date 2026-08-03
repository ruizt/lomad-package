# Internal fit implementation for lomad_fit()
#
# .lomad_fit_clt() — CLT pipeline (paper method)


# CLT pipeline: the paper method.
# Estimates all quantities needed for pointwise inference on R_t.
.lomad_fit_clt <- function(y1, y2, h, s, lag_max,
                          noise_method = "ar1", noise_override = NULL,
                          min_lambda = 0) {

  y1 <- as.numeric(y1)
  y2 <- as.numeric(y2)
  n  <- length(y1)

  if (is.null(h)) h <- max(5L, floor(n / 200L))
  if (is.null(s)) {
    s <- min(60L * as.integer(h), floor(n / 4L))
    message(sprintf("h = %d, s = %d (auto)", h, s))
  }
  h       <- as.integer(h)
  s       <- as.integer(s)
  lag_max <- as.integer(lag_max)
  if (s <= 3L) stop("`s` must be > 3.")
  if (h < 1L)  stop("`h` must be >= 1.")

  tr <- estimate_trends(y1, y2, h)

  raw_lag_max <- lag_max + h - 1L

  if (!is.null(noise_override)) {
    # Accept a single spec (shared) or a list of two (per-series).
    # A spec is either parametric (ar, ma, sigma2) or a raw ACVF (acov).
    is_spec <- function(x)
      is.list(x) && (!is.null(x$ar) || !is.null(x$acov))
    if (is_spec(noise_override)) {
      spec1 <- spec2 <- noise_override
    } else {
      spec1 <- noise_override[[1]]
      spec2 <- noise_override[[2]]
    }
    norm_spec <- function(sp) {
      if (!is.null(sp$acov)) return(list(acov = as.numeric(sp$acov)))
      list(ar = sp$ar,
           ma = if (!is.null(sp$ma)) sp$ma else numeric(0),
           sigma2 = sp$sigma2)
    }
    noise <- list(series1 = norm_spec(spec1), series2 = norm_spec(spec2))
    message("Using noise_override (oracle parameters)")
  } else {
    noise <- estimate_ar1_noise(y1, y2, tr$ma1, tr$ma2, h = h)
    message(sprintf("AR(1): phi = %.3f / %.3f",
                    noise$series1$ar, noise$series2$ar))
  }

  # Raw noise ACVF to lag raw_lag_max: from the parametric model, or directly
  # from a nonparametric/override ACVF (zero-padded past its length).
  raw_acov <- function(sp) {
    if (!is.null(sp$acov)) {
      out <- numeric(raw_lag_max + 1L)
      k   <- min(length(sp$acov), raw_lag_max + 1L)
      out[seq_len(k)] <- sp$acov[seq_len(k)]
      return(out)
    }
    arma_acov(ar = sp$ar, ma = sp$ma, sigma2 = sp$sigma2,
              lag_max = raw_lag_max)
  }
  acov1_raw <- raw_acov(noise$series1)
  acov2_raw <- raw_acov(noise$series2)

  acov1 <- .ma_filter_acov(acov1_raw, h, lag_max)
  acov2 <- .ma_filter_acov(acov2_raw, h, lag_max)
  sums  <- acov_sums(acov1, acov2)

  # Long-run variances are nonnegative for any valid spectral density;
  # nonparametric estimates can dip below zero from sampling noise.
  for (nm in c("L1", "L2")) {
    if (sums[[nm]] < 0) {
      warning(sprintf(
        "Estimated long-run sum %s was negative (%.3g); floored at 0.",
        nm, sums[[nm]]))
      sums[[nm]] <- 0
    }
  }

  sigma1_sq <- acov1[1L]
  sigma2_sq <- acov2[1L]

  # Each smoothed series carries its own MA-filtered noise, so the windowed
  # sample variance of ma_k over-states tau_k^2 by the *expected windowed*
  # variance of that noise -- not its marginal variance, which over-subtracts
  # when the smoothed noise is autocorrelated (the window mean absorbs
  # low-frequency variation). Under the old shared-trend scheme the averaging
  # in (ma1 + ma2)/2 halved the noise before tau^2 was formed; per-series
  # estimation loses that, so this correction is load-bearing rather than
  # cosmetic.
  bias1 <- max(0, .windowed_var_expect(acov1, s))
  bias2 <- max(0, .windowed_var_expect(acov2, s))

  tau1_sq <- pmax(0, compute_tau_sq(tr$ma1, s) - bias1)
  tau2_sq <- pmax(0, compute_tau_sq(tr$ma2, s) - bias2)

  rho <- compute_rho(tau1_sq, tau2_sq, sigma1_sq, sigma2_sq)   # rho^(0), r = 1
  V   <- compute_V(tau1_sq, tau2_sq, sigma1_sq, sigma2_sq,
                   sums$L1, sums$L2, sums$Q1, sums$Q2, sums$Q12)

  R <- rep(NA_real_, n)
  for (t in s:n) {
    w   <- (t - s + 1L):t
    y1w <- tr$ma1[w]
    y2w <- tr$ma2[w]
    if (any(is.na(y1w)) || any(is.na(y2w))) next
    if (stats::sd(y1w) == 0 || stats::sd(y2w) == 0) next
    r <- suppressWarnings(stats::cor(y1w, y2w))
    if (is.finite(r) && abs(r) < 1) R[t] <- r
  }

  # Per-window signal-to-noise. rho^(0) = tau1*tau2/sqrt(BD) goes to zero as
  # either lambda does, so where one series carries almost no trend signal the
  # null predicts near-zero correlation and the test has no power: R has
  # nothing to fall below. Those windows are *untestable*, which is a
  # different statement from "no decoupling detected", and reporting them as
  # the latter is misleading. Note this is an information limit, not an
  # artefact of estimating two taus -- at lambda near zero a small trend
  # proportional to the other series and a flat trend produce the same
  # observed correlation.
  lambda1  <- tau1_sq / sigma1_sq
  lambda2  <- tau2_sq / sigma2_sq
  testable <- is.finite(lambda1) & is.finite(lambda2) &
    pmin(lambda1, lambda2) >= min_lambda

  # r_hat = R / rho^(0) is the affine-invariant effect-size estimate: it is
  # the sample correlation with the attenuation from finite SNR divided out,
  # so delta_hat = sqrt(1 - r_hat^2) estimates the local separation directly.
  r_hat <- ifelse(is.finite(rho) & rho > 0, R / rho, NA_real_)

  valid_idx <- which(is.finite(R) & is.finite(rho) & is.finite(V) & V > 0 &
                     testable)

  list(
    method    = "clt",
    ma1       = tr$ma1,
    ma2       = tr$ma2,
    noise     = noise,
    acov_sums = sums,
    tau1_sq   = tau1_sq,
    tau2_sq   = tau2_sq,
    rho       = rho,
    V         = V,
    R         = R,
    r_hat     = r_hat,
    lambda1   = lambda1,
    lambda2   = lambda2,
    testable  = testable,
    valid_idx = valid_idx,
    inputs    = list(n = n, h = h, s = s, lag_max = lag_max)
  )
}
