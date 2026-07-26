#' Estimate AR(1) noise parameters via the variogram
#'
#' Estimates the AR(1) coefficient and innovation variance for each series
#' using the lag-1 and lag-2 sample variograms of shared-trend residuals.
#' The variogram at short lags is dominated by the noise (not the trend),
#' making the estimates robust to trend differences when \eqn{d > 0}.
#'
#' This is a special case of the difference-based approach of
#' Hall and Van Keilegom (2003) for the AR(1) model. For general AR(p)
#' noise, see [estimate_arma_noise()].
#'
#' Residuals are computed as \eqn{y_k - \hat s_t}, where \eqn{\hat s_t} is
#' the shared trend estimate from [estimate_trends()].  The variogram
#' \deqn{V(l) = \frac{1}{2(n - l)} \sum_{t=1}^{n-l} (r_{t+l} - r_t)^2}
#' satisfies \eqn{V(l) = \gamma(0) - \gamma(l)} for a stationary process.
#' For AR(1) with coefficient \eqn{\phi} and innovation variance
#' \eqn{\sigma^2}: \eqn{V(1) = \sigma^2 / (1 + \phi)} and
#' \eqn{V(2) = \sigma^2}, giving \eqn{\phi = V(2)/V(1) - 1}.
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param trend Numeric vector (length `n`). Shared trend estimate
#'   \eqn{\hat s_t} as returned by [estimate_trends()]`$trend`.
#'
#' @return A named list with elements `series1` and `series2`, each containing:
#'   \describe{
#'     \item{ar}{Numeric scalar. Estimated AR(1) coefficient \eqn{\hat\phi}.}
#'     \item{ma}{`numeric(0)` (no MA component).}
#'     \item{sigma2}{Positive numeric. Innovation variance \eqn{\hat\sigma^2}.}
#'     \item{order}{Integer vector `c(1, 0, 0)`.}
#'   }
#'
#' @references
#' Hall, P. and Van Keilegom, I. (2003). Using difference-based methods for
#' inference in nonparametric regression with time series errors. \emph{Journal
#' of the Royal Statistical Society Series B}, 65(2), 443--456.
#'
#' @seealso [estimate_arma_noise()] for general AR(p) noise estimation.
#'
#' @export
estimate_ar1_noise <- function(y1, y2, trend) {
  if (length(y1) != length(y2) || length(y1) != length(trend))
    stop("`y1`, `y2`, and `trend` must all have the same length.")

  valid_t <- which(!is.na(trend))
  if (length(valid_t) < 10L)
    stop("Fewer than 10 valid trend points; cannot estimate noise.")

  resid1 <- y1[valid_t] - trend[valid_t]
  resid2 <- y2[valid_t] - trend[valid_t]

  list(
    series1 = .variogram_ar1(resid1),
    series2 = .variogram_ar1(resid2)
  )
}


#' Estimate AR(p) noise parameters via the multi-lag variogram
#'
#' Generalises [estimate_ar1_noise()] to AR(p) noise using multi-lag
#' variograms and Yule-Walker estimation, following the difference-based
#' framework of Hall and Van Keilegom (2003).
#'
#' The variogram \eqn{V(l) = \gamma(0) - \gamma(l)} is computed at lags
#' \eqn{1, \ldots, L}.  At sufficiently large \eqn{L}, \eqn{V(L) \approx
#' \gamma(0)}, giving an estimate of the process variance.  The recovered
#' autocovariances \eqn{\hat\gamma(l) = \hat\gamma(0) - V(l)} are then
#' plugged into the Yule-Walker equations to obtain AR coefficients and the
#' innovation variance.
#'
#' If `p` is supplied, the AR order is fixed.
#' If `p = NULL` (default), BIC selects among orders \eqn{1, \ldots, p_{max}}.
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param trend Numeric vector (length `n`). Shared trend estimate
#'   \eqn{\hat s_t} as returned by [estimate_trends()]`$trend`.
#' @param p Integer or `NULL`. If specified, fit AR(p) directly (no order
#'   selection). If `NULL`, select order via BIC.
#' @param p_max Positive integer. Maximum AR order for BIC selection.
#'   Ignored when `p` is specified. Default 5.
#' @param L Positive integer. Maximum variogram lag used to estimate
#'   \eqn{\gamma(0)}. Default 20.
#'
#' @return A named list with elements `series1` and `series2`, each containing:
#'   \describe{
#'     \item{ar}{Numeric vector of AR coefficients.}
#'     \item{ma}{`numeric(0)` (no MA component in the AR approximation).}
#'     \item{sigma2}{Positive numeric. Innovation variance \eqn{\hat\sigma^2}.}
#'     \item{order}{Integer vector `c(p, 0, 0)`.}
#'   }
#'
#' @references
#' Hall, P. and Van Keilegom, I. (2003). Using difference-based methods for
#' inference in nonparametric regression with time series errors. \emph{Journal
#' of the Royal Statistical Society Series B}, 65(2), 443--456.
#'
#' @seealso [estimate_ar1_noise()] for the AR(1) special case.
#'
#' @export
estimate_arma_noise <- function(y1, y2, trend,
                                p     = NULL,
                                p_max = 5L,
                                L     = 20L) {
  if (length(y1) != length(y2) || length(y1) != length(trend))
    stop("`y1`, `y2`, and `trend` must all have the same length.")

  p_max <- as.integer(p_max)
  L     <- as.integer(L)
  if (!is.null(p)) p <- as.integer(p)

  valid_t <- which(!is.na(trend))
  if (length(valid_t) < 10L)
    stop("Fewer than 10 valid trend points; cannot estimate noise.")

  resid1 <- y1[valid_t] - trend[valid_t]
  resid2 <- y2[valid_t] - trend[valid_t]

  list(
    series1 = .variogram_ar(resid1, p = p, p_max = p_max, L = L),
    series2 = .variogram_ar(resid2, p = p, p_max = p_max, L = L)
  )
}


#' Estimate noise autocovariance nonparametrically via the variogram
#'
#' Model-free alternative to [estimate_ar1_noise()] and
#' [estimate_arma_noise()]: recovers the noise autocovariance function
#' directly from the multi-lag variogram of shared-trend residuals, without
#' fitting an AR or ARMA model. Use this as a robustness fallback when the
#' parametric noise models are misspecified (e.g. oscillatory or
#' near-noninvertible noise) — situations the parametric estimators signal
#' with boundary warnings or implausible variance estimates.
#'
#' The lag-\eqn{l} variogram \eqn{V(l)} of the residuals satisfies
#' \eqn{V(l) = \gamma(0) - \gamma(l)} for stationary noise and is robust to
#' trend contamination (Hall and Van Keilegom, 2003). The recovered sequence
#' \eqn{\hat\gamma(l) = \hat\gamma(0) - V(l)} is Bartlett-tapered (weights
#' \eqn{(1 - l/(b+1))_+}) so that the long-run sums entering the asymptotic
#' variance are stable, in the spirit of Newey and West (1987).
#'
#' The variance anchor \eqn{\hat\gamma(0)} is controlled by `anchor`. The
#' default `"variance"` uses the residual sample variance: a tight estimator
#' that keeps the downstream one-sided test calibrated; residual trend
#' contamination can only inflate it, biasing the test in the conservative
#' direction. `"plateau"` uses the variogram plateau average: fully
#' trend-robust, but its sampling error propagates coherently to every lag
#' of \eqn{\hat\gamma}, which can make the test anti-conservative — prefer
#' it only for estimation (not testing) under strong suspected separation.
#'
#' Compared with the parametric estimators, this method removes model-class
#' risk at the cost of higher variance in the recovered autocovariances
#' (especially at long lags) and one tuning constant (the taper bandwidth).
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param trend Numeric vector (length `n`). Shared trend estimate
#'   \eqn{\hat s_t} as returned by [estimate_trends()]`$trend`.
#' @param lag_max Positive integer. Length of the returned autocovariance
#'   sequence (lags 0 to `lag_max`; zero beyond the taper bandwidth).
#'   Default 100.
#' @param bandwidth Positive integer or `NULL`. Bartlett taper bandwidth
#'   \eqn{b}. `NULL` (default) uses \eqn{\lceil 10 \log_{10} n \rceil}.
#' @param anchor Character. \eqn{\gamma(0)} anchor: `"variance"` (default,
#'   residual sample variance — calibrated testing, conservative under
#'   separation) or `"plateau"` (variogram plateau — trend-robust but can
#'   make testing anti-conservative). See Details.
#'
#' @return A named list with elements `series1` and `series2`, each
#'   containing:
#'   \describe{
#'     \item{acov}{Numeric vector of length `lag_max + 1`. Estimated raw
#'       noise autocovariances \eqn{\hat\gamma(0), \ldots,
#'       \hat\gamma(\text{lag\_max})}.}
#'   }
#'
#' @references
#' Hall, P. and Van Keilegom, I. (2003). Using difference-based methods for
#' inference in nonparametric regression with time series errors. \emph{Journal
#' of the Royal Statistical Society Series B}, 65(2), 443--456.
#'
#' Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite,
#' heteroskedasticity and autocorrelation consistent covariance matrix.
#' \emph{Econometrica}, 55(3), 703--708.
#'
#' @seealso [estimate_ar1_noise()], [estimate_arma_noise()], [lomad_fit()]
#'
#' @examples
#' cpl <- subset(sim_decoupling, scenario == "coupled")
#' tr  <- estimate_trends(cpl$y1, cpl$y2, h = 5)
#' nz  <- estimate_acf_noise(cpl$y1, cpl$y2, tr$trend, lag_max = 30)
#' plot(0:30, nz$series1$acov, type = "h", xlab = "lag",
#'      ylab = "estimated noise ACVF")
#'
#' @export
estimate_acf_noise <- function(y1, y2, trend,
                               lag_max   = 100L,
                               bandwidth = NULL,
                               anchor    = c("variance", "plateau")) {
  if (length(y1) != length(y2) || length(y1) != length(trend))
    stop("`y1`, `y2`, and `trend` must all have the same length.")

  lag_max <- as.integer(lag_max)
  anchor  <- match.arg(anchor)

  valid_t <- which(!is.na(trend))
  if (length(valid_t) < 30L)
    stop("Fewer than 30 valid trend points; cannot estimate the noise ACVF.")

  resid1 <- y1[valid_t] - trend[valid_t]
  resid2 <- y2[valid_t] - trend[valid_t]

  list(
    series1 = list(acov = .variogram_acov(resid1, lag_max, bandwidth, anchor)),
    series2 = list(acov = .variogram_acov(resid2, lag_max, bandwidth, anchor))
  )
}
