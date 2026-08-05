# CLT building blocks: autocovariance filtering and asymptotic quantities
#
# Internal:
#   .ma_filter_acov()  — ACVF of the MA(h)-smoothed process
#
# Exported (used by simulation scripts and lomad_fit_clt pipeline):
#   compute_tau_sq()   — rolling signal variance
#   compute_rho()      — local population correlation
#   compute_V()        — asymptotic variance of local sample correlation


# ---- Internal helpers ----------------------------------------------------

# ACVF of the MA(h)-smoothed process.
#
# If Z_t has ACVF gamma_Z(l), then eta_t = (1/h) * sum_{j=0}^{h-1} Z_{t-j}
# has ACVF:
#   gamma_eta(l) = (1/h^2) * sum_{m=-(h-1)}^{h-1} (h - |m|) * gamma_Z(|l + m|)
#
# acov_raw : gamma_Z(0), gamma_Z(1), ... of length >= lag_max + h
# h        : MA smoothing window width
# lag_max  : maximum lag to compute for eta
.ma_filter_acov <- function(acov_raw, h, lag_max) {
  m_vals <- -(h - 1L):(h - 1L)
  wts    <- (h - abs(m_vals)) / h^2
  result <- numeric(lag_max + 1L)
  for (l in 0L:lag_max) {
    lag_vals       <- abs(l + m_vals)
    result[l + 1L] <- sum(wts * acov_raw[lag_vals + 1L])
  }
  result
}

# Expected windowed sample variance of a stationary process.
#
# For a stationary process xi_t with ACVF gamma and a window W of s points,
# the population-denominator sample variance satisfies
#   E[(1/s) sum_{t in W} (xi_t - xi_bar_W)^2]
#     = gamma(0) - (1/s^2) sum_{t,u in W} gamma(t - u)
#     = (1 - 1/s) gamma(0) - (2/s^2) sum_{l=1}^{s-1} (s - l) gamma(l)
#
# This is smaller than the marginal variance gamma(0) whenever the process is
# positively autocorrelated, because the window mean absorbs low-frequency
# variation. Lags beyond length(acov) - 1 are treated as zero (truncation).
#
# acov : gamma(0), gamma(1), ... (nonnegative lags)
# s    : window length
.windowed_var_expect <- function(acov, s) {
  s <- as.integer(s)
  stopifnot(s >= 2L, length(acov) >= 1L)
  l <- seq_len(min(s - 1L, length(acov) - 1L))
  (1 - 1 / s) * acov[1L] - (2 / s^2) * sum((s - l) * acov[l + 1L])
}


# ---- Exported functions --------------------------------------------------

#' Rolling signal variance (tau-squared) in each window
#'
#' For each time point \eqn{t}, computes the variance of the trend over the
#' left-aligned window \eqn{W_t = \{t - s + 1, \ldots, t\}}:
#' \deqn{\hat\tau^2_t = \frac{1}{s} \sum_{u \in W_t} (s_u - \bar s)^2.}
#' Returns `NA` for positions where fewer than `s` trend values are available.
#'
#' @param trend Numeric vector. Estimated (or true) shared trend.
#' @param s Positive integer. Rolling window length (paper notation: \eqn{s}).
#'
#' @return Numeric vector of length `length(trend)` with `NA` for the first
#'   `s - 1` positions.
#'
#' @export
compute_tau_sq <- function(trend, s) {
  s <- as.integer(s)
  stopifnot(s >= 2L, is.numeric(trend))
  n      <- length(trend)
  tau_sq <- rep(NA_real_, n)
  for (t in s:n) {
    w <- trend[(t - s + 1L):t]
    if (any(is.na(w))) next
    tau_sq[t] <- mean((w - mean(w))^2)
  }
  tau_sq
}


#' Local population correlation under affine similarity
#'
#' Computes the approximation for the local population correlation \eqn{\rho}
#' from Proposition 1:
#' \deqn{\rho \approx
#'   \frac{r\,\tau_1\tau_2}{\sqrt{(\tau_1^2 + \sigma_1^2)(\tau_2^2 + \sigma_2^2)}}.}
#' The default \eqn{r = 1} gives \eqn{\rho^{(0)}}, the value under \eqn{H_0},
#' which is what the test uses. Values \eqn{r < 1} give the attenuated correlation
#' under the alternative, via the exact identity
#' \eqn{\rho = (1 - \delta^2)^{1/2}\rho^{(0)}} with \eqn{\delta = (1-r^2)^{1/2}}.
#'
#' @param tau1_sq Non-negative numeric (scalar or vector). Signal variance
#'   \eqn{\tau_1^2} of series 1 over the local window.
#' @param tau2_sq Non-negative numeric (scalar or vector). Signal variance
#'   \eqn{\tau_2^2} of series 2 over the local window.
#' @param sigma1_sq Positive numeric (scalar). Marginal noise variance for
#'   series 1, \eqn{\sigma_1^2 = \gamma_1(0)}.
#' @param sigma2_sq Positive numeric (scalar). Marginal noise variance for
#'   series 2, \eqn{\sigma_2^2 = \gamma_2(0)}.
#' @param r Numeric in \eqn{[-1, 1]}. Window correlation of the two signals.
#'   Defaults to 1 (the null).
#'
#' @return Numeric (recycled to the length of `tau1_sq`/`tau2_sq`).
#'
#' @export
compute_rho <- function(tau1_sq, tau2_sq, sigma1_sq, sigma2_sq, r = 1) {
  stopifnot(sigma1_sq > 0, sigma2_sq > 0,
            all(tau1_sq >= 0, na.rm = TRUE),
            all(tau2_sq >= 0, na.rm = TRUE))
  B <- tau1_sq + sigma1_sq
  D <- tau2_sq + sigma2_sq
  r * sqrt(tau1_sq * tau2_sq) / sqrt(B * D)
}


#' Asymptotic variance of the local sample correlation (V)
#'
#' Computes the approximation for \eqn{V} from Proposition 1 under \eqn{H_0}.
#' The pointwise CLT gives \eqn{\sqrt{s}(R_t - \rho_t) \to N(0, V_t)}, so the
#' standard error of \eqn{R_t} is \eqn{\sqrt{V_t / s}}.
#'
#' Note the cross-pairing: \eqn{\tau_2^2} multiplies the \eqn{L_1} term and
#' \eqn{\tau_1^2} multiplies the \eqn{L_2} term. This comes from the
#' \eqn{(5,5)} entry of \eqn{\Sigma}, where
#' \eqn{\mathrm{Cov}(\tilde Y_{1t}\tilde Y_{2t}, \cdot)} contributes
#' \eqn{s_{1t}s_{1,t-l}\gamma_2(l) + s_{2t}s_{2,t-l}\gamma_1(l)}.
#'
#' @param tau1_sq Non-negative numeric (scalar or vector). Signal variance,
#'   series 1.
#' @param tau2_sq Non-negative numeric (scalar or vector). Signal variance,
#'   series 2.
#' @param sigma1_sq Positive numeric. Marginal noise variance, series 1.
#' @param sigma2_sq Positive numeric. Marginal noise variance, series 2.
#' @param L1 Numeric. \eqn{L_1 = \sum_{l} \gamma_1(l)}, long-run variance of
#'   series 1 noise.
#' @param L2 Numeric. \eqn{L_2 = \sum_{l} \gamma_2(l)}.
#' @param Q1 Numeric. \eqn{Q_1 = \sum_{l} \gamma_1(l)^2}.
#' @param Q2 Numeric. \eqn{Q_2 = \sum_{l} \gamma_2(l)^2}.
#' @param Q12 Numeric. \eqn{Q_{12} = \sum_{l} \gamma_1(l)\gamma_2(l)}.
#'
#' @return Numeric, non-negative.
#'
#' @export
compute_V <- function(tau1_sq, tau2_sq, sigma1_sq, sigma2_sq,
                      L1, L2, Q1, Q2, Q12) {
  stopifnot(sigma1_sq > 0, sigma2_sq > 0)
  B <- tau1_sq + sigma1_sq
  D <- tau2_sq + sigma2_sq
  Q12 / (B * D) +
    tau2_sq * sigma1_sq^2 * L1 / (B^3 * D) +
    tau1_sq * sigma2_sq^2 * L2 / (B * D^3) +
    tau1_sq * tau2_sq * Q1 / (2 * B^3 * D) +
    tau1_sq * tau2_sq * Q2 / (2 * B * D^3)
}
