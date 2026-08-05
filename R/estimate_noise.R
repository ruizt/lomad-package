#' Estimate AR(1) noise parameters via the variogram
#'
#' Estimates the AR(1) coefficient and innovation variance for each series
#' from the lag-1 and lag-2 sample variograms. Differencing annihilates a
#' slowly varying trend, so the variogram at short lags is dominated by the
#' noise and no trend estimate is required.
#'
#' This is a special case of the difference-based approach of
#' Hall and Van Keilegom (2003) for the AR(1) model. For general AR(p)
#'
#' Each series is detrended by its own smoothed version before differencing.
#' Differencing alone annihilates a linear trend but not a curved one, and at
#' the smoothing windows used in practice the residual curvature inflates the
#' short-lag variogram enough to bias \eqn{\hat\phi} upward; subtracting a
#' shared trend instead makes that bias grow with the separation between the
#' two trends. Per-series detrending avoids both. The variogram
#' \deqn{V(l) = \frac{1}{2(n - l)} \sum_{t=1}^{n-l} (r_{t+l} - r_t)^2}
#' satisfies \eqn{V(l) = \gamma(0) - \gamma(l)} for a stationary process.
#' For AR(1) with coefficient \eqn{\phi} and innovation variance
#' \eqn{\sigma^2}: \eqn{V(1) = \sigma^2 / (1 + \phi)} and
#' \eqn{V(2) = \sigma^2}, giving \eqn{\phi = V(2)/V(1) - 1}.
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param trend1 Optional numeric vector (length `n`). Trend for series 1,
#'   subtracted before differencing. `NULL` differences `y1` directly.
#' @param trend2 Optional numeric vector (length `n`). As `trend1`, for
#'   series 2. Passing the same vector for both reproduces the pre-0.1.0
#'   shared-trend behaviour.
#' @param h Optional integer. Width of the moving average subtracted to form
#'   the residuals. Supplying it corrects the estimate for the high-pass
#'   filtering that subtraction induces; `NULL` skips the correction and is
#'   appropriate only when no trend was subtracted.
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
estimate_ar1_noise <- function(y1, y2, trend1 = NULL, trend2 = NULL,
                               h = NULL) {
  if (length(y1) != length(y2))
    stop("`y1` and `y2` must have the same length.")
  for (nm in c("trend1", "trend2")) {
    tr <- get(nm)
    if (!is.null(tr) && length(tr) != length(y1))
      stop(sprintf("`%s` must have the same length as `y1` and `y2`.", nm))
  }

  r1 <- if (is.null(trend1)) y1 else y1 - trend1
  r2 <- if (is.null(trend2)) y2 else y2 - trend2

  valid_t <- which(is.finite(r1) & is.finite(r2))
  if (length(valid_t) < 10L)
    stop("Fewer than 10 valid residual points; cannot estimate noise.")
  resid1 <- r1[valid_t]
  resid2 <- r2[valid_t]

  list(
    series1 = .variogram_ar1(resid1, h = h),
    series2 = .variogram_ar1(resid2, h = h)
  )
}
