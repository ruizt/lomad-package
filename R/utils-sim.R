# Internal utilities for sim_trends() and sim_noise()
#
# Fourier basis:
#   .generate_fourier_coef()  — draw coefficients with spectral decay
#   .generate_coef_pair()     — pair of coefficient vectors at distance d
#   .make_basis_trends()      — evaluate Fourier basis to get mu1, mu2, x_mean
#
# Coupling weights:
#   .make_w_smooth()          — stochastic repulsion
#   .make_w_cross()           — stochastic crossing
#   .make_w_rate()            — periodic event-based decoupling
#
# Mixing:
#   .apply_w()                — mix mu1, mu2 via w and rescale to target d
#
# Noise:
#   .pacf_to_arma_coefs()     — PACF parameterisation → ARMA coefficients


# ---- Fourier basis -------------------------------------------------------

# Draw Fourier coefficients with spectral decay: sd_k = sd0 / k^p
.generate_fourier_coef <- function(nb, sd0 = 2, p = 2.5, k_min = 1L,
                                   normalize = TRUE) {
  K    <- (nb - 1L) / 2L
  k    <- rep(seq_len(K), each = 2L)
  sd_k <- ifelse(k >= k_min, sd0 / (k^p), 0)
  cf   <- stats::rnorm(2L * K, mean = 0, sd = sd_k)

  # Fix the total signal power, keeping the spectral shape random.
  #
  # The affine effect size is delta = d / sqrt(||mu1||^2 + d^2), so a random
  # ||mu1|| feeds straight into it: unnormalized, ||mu1|| has CV 0.51 and
  # spans 0.34 to 7.41 across seeds, which spreads delta from roughly 0.2 to
  # 0.7 at a nominal d = 1. Each point on a power curve would then average
  # over a wide band of true effect sizes.
  #
  # Amplitude is not a design factor worth keeping random here: SNR is already
  # controlled separately via lambda_target, and the trend *shape* still
  # varies freely. Normalizing to sqrt(sum(sd_k^2)) = E||coef||^2 ^ (1/2)
  # leaves the expected scale (and hence the meaning of sd0) untouched and
  # removes only the fluctuation.
  if (normalize) {
    target <- sqrt(sum(sd_k^2))
    nrm    <- sqrt(sum(cf^2))
    if (nrm > 0 && target > 0) cf <- cf * (target / nrm)
  }
  cf
}

# Generate a pair of coefficient vectors separated by distance d using the
# sphere-to-ellipse transform method.
.generate_coef_pair <- function(nb, sd0 = 2, d = 1, p = 2.5, k_min = 1L,
                                seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  coef1 <- .generate_fourier_coef(nb, sd0, p, k_min)

  z <- stats::rnorm(length(coef1))
  u <- z / sqrt(sum(z^2))

  k_idx <- rep(seq_len((nb - 1L) / 2L), each = 2L)
  axes  <- ifelse(k_idx >= k_min, 1 / (k_idx^max(0, p - 0.1 * d)), 0)
  dir   <- u * axes

  # Project out the component along coef1. Under affine similarity H_0 holds
  # iff coef2 is parallel to coef1, so the null is a ray, not a point, and a
  # displacement along coef1 moves coef2 *within* the null set -- it raises d
  # without raising the departure from H_0. Leaving it in makes d a
  # mismeasurement of effect size rather than a noisy measure of it.
  c1_norm2 <- sum(coef1^2)
  if (c1_norm2 > 0) dir <- dir - (sum(dir * coef1) / c1_norm2) * coef1

  dn <- sqrt(sum(dir^2))
  if (dn < .Machine$double.eps^0.5)
    stop("Displacement direction collapsed after projection; redraw with a ",
         "different seed.")
  dir <- dir / dn

  coef2 <- coef1 + d * dir

  # Rescale coef2 to ||coef1||. The displacement is orthogonal, so without
  # this ||mu2||^2 = ||mu1||^2 + d^2 and the two series carry different signal
  # power -- lambda_2/lambda_1 reaches 3 by d = 3, drifting with the very axis
  # the power curves are plotted against.
  #
  # This costs nothing. The effect size is delta = sqrt(1 - r^2) with
  # r = cos(angle between the trends), and rescaling is a pure scale change,
  # which is exactly what the affine null is invariant to. So r, and hence
  # delta = d/sqrt(||mu1||^2 + d^2), are preserved exactly. What it gives up
  # is ||mu1 - mu2|| = d, which under a scale-invariant null is no longer the
  # quantity worth preserving.
  c2_norm2 <- sum(coef2^2)
  if (c2_norm2 > 0) coef2 <- coef2 * sqrt(c1_norm2 / c2_norm2)

  list(coef1 = coef1, coef2 = coef2)
}

# Evaluate Fourier basis pair: returns mu1, mu2, x_mean
.make_basis_trends <- function(n, nb, sd0, d, p, k_min, seed) {
  coefs <- .generate_coef_pair(nb = nb, sd0 = sd0, d = d, p = p,
                               k_min = k_min, seed = seed)

  fb  <- fda::create.fourier.basis(rangeval = c(0, n), nbasis = nb, period = n)
  Phi <- fda::eval.basis(seq_len(n), fb)[, -1]

  mu1    <- as.numeric(Phi %*% coefs$coef1)
  mu2    <- as.numeric(Phi %*% coefs$coef2)
  x_mean <- (mu1 + mu2) / 2

  list(mu1 = mu1, mu2 = mu2, x_mean = x_mean)
}


# ---- Coupling weight generators -----------------------------------------

# Stochastic repulsion: w in (0, 1), no crossing.
# coupling controls the fraction of time w > 0.5.
.make_w_smooth <- function(n, bw = 50, coupling = 0.8) {
  if (coupling <= 0 || coupling >= 1) stop("`coupling` must be in (0, 1).")
  z_raw    <- stats::rnorm(n)
  z_smooth <- stats::ksmooth(seq_len(n), z_raw, kernel = "normal",
                             bandwidth = bw, x.points = seq_len(n))$y
  z <- (z_smooth - mean(z_smooth)) / stats::sd(z_smooth)
  stats::pnorm(stats::qnorm(coupling) - z)
}

# Stochastic crossing: w centered at 1, can go below 0 or above 1.
# coupling controls the fraction of time |w - 1| < 0.5.
.make_w_cross <- function(n, bw = 50, coupling = 0.8) {
  if (coupling <= 0 || coupling >= 1) stop("`coupling` must be in (0, 1).")
  z_raw    <- stats::rnorm(n)
  z_smooth <- stats::ksmooth(seq_len(n), z_raw, kernel = "normal",
                             bandwidth = bw, x.points = seq_len(n))$y
  z        <- (z_smooth - mean(z_smooth)) / stats::sd(z_smooth)
  sigma_w  <- 0.5 / stats::qnorm((1 + coupling) / 2)
  1 - sigma_w * z
}

# Periodic event-based decoupling: evenly spaced dips at fixed times.
#
# `bump` sets the pulse shape and changes nothing else -- event times, spacing
# and depth are identical either way.
#
#   "gaussian" (default) symmetric normal. Smooth everywhere; no onset corner.
#   "gamma"    shape-2 gamma, matched on width. Rises from zero with non-zero
#              slope, so w_t has a corner at each event onset.
#
# Gaussian is the default because a corner is exactly what this package's own
# noise estimation cannot handle: it is difference-based, and trend-robust only
# for trends with bounded derivative (Hall and Van Keilegom, 2003, eqn 2.4). The
# gamma pulse therefore leaks into the residual autocovariance far more than the
# gaussian one, badly enough to cost the fixed-rate structure most of its
# detection at high autocorrelation. Gamma is kept because that contrast is
# itself worth simulating, not as a sensible starting point.
.make_w_rate <- function(n, rate = 0.01, bump = c("gaussian", "gamma")) {
  bump     <- match.arg(bump)
  n_events <- round(rate * n)
  if (n_events < 1) stop("`rate * n` must be at least 1; increase `rate` or `n`.")

  delta     <- 1 / n_events
  gap       <- n / n_events
  dl        <- round(gap * delta)
  strength  <- n_events

  scale_g <- dl * 0.5
  sd_g    <- scale_g * sqrt(2)   # sd of gamma(shape = 2, scale = scale_g)

  events <- round(seq(from = gap / 2, by = gap, length.out = n_events))
  t      <- seq_len(n)
  w      <- rep(1, n)
  for (i in events) {
    b <- switch(bump,
      gamma    = stats::dgamma(t - i, shape = 2, scale = scale_g),
      gaussian = stats::dnorm(t - i - sd_g, sd = sd_g)
    )
    b <- b / max(b)
    w <- w - delta * strength * b
  }
  pmax(0, pmin(1, w))
}


# ---- Mixing and rescaling -----------------------------------------------

# Mix mu1, mu2 via coupling weight w and rescale to target distance d.
#
# Mixing formula: x_i = x_mean + (1 - w) * (mu_i - x_mean)
# Then rescale so that ||x1 - x2|| = d.
.apply_w <- function(mu1, mu2, x_mean, w, d) {
  x1_unit <- x_mean + (1 - w) * (mu1 - x_mean)
  x2_unit <- x_mean + (1 - w) * (mu2 - x_mean)

  r_unit <- sqrt(sum((x1_unit - x2_unit)^2))
  s      <- d / r_unit

  list(
    x1     = x_mean + s * (x1_unit - x_mean),
    x2     = x_mean + s * (x2_unit - x_mean),
    x_mean = x_mean,
    w      = w
  )
}


# ---- Noise utilities -----------------------------------------------------

# Convert p PACF values (drawn uniformly from (-0.95, 0.95)) to ARMA
# coefficients via the Durbin-Levinson recursion. All resulting processes
# are causal/invertible by construction.
.pacf_to_arma_coefs <- function(p) {
  if (p == 0) return(numeric(0))
  pacf <- stats::runif(p, min = -0.95, max = 0.95)
  coefs <- pacf[1]
  if (p == 1) return(coefs)
  for (k in 2:p) {
    coefs <- c(coefs - pacf[k] * rev(coefs), pacf[k])
  }
  coefs
}
