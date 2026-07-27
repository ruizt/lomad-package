#' Estimate AR(1) noise parameters via the variogram
#'
#' Estimates the AR(1) coefficient and innovation variance for each series
#' using the lag-1 and lag-2 sample variograms of shared-trend residuals.
#' The variogram at short lags is dominated by the noise (not the trend),
#' making the estimates robust to trend differences when \eqn{d > 0}.
#'
#' This is a special case of the difference-based approach of
#' Hall and Van Keilegom (2003) for the AR(1) model. For general AR(p)
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
