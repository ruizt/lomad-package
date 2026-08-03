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
#' By default the variogram is taken on the observed series directly. Passing
#' `trend` instead subtracts it first, reproducing the pre-0.1.0 behaviour;
#' that path is retained only for comparison, since under affine similarity
#' there is no shared trend to subtract. The variogram
#' \deqn{V(l) = \frac{1}{2(n - l)} \sum_{t=1}^{n-l} (r_{t+l} - r_t)^2}
#' satisfies \eqn{V(l) = \gamma(0) - \gamma(l)} for a stationary process.
#' For AR(1) with coefficient \eqn{\phi} and innovation variance
#' \eqn{\sigma^2}: \eqn{V(1) = \sigma^2 / (1 + \phi)} and
#' \eqn{V(2) = \sigma^2}, giving \eqn{\phi = V(2)/V(1) - 1}.
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param trend Optional numeric vector (length `n`). If supplied, residuals
#'   \eqn{y_k - \hat s_t} are formed before differencing. Defaults to `NULL`,
#'   which differences the observed series directly.
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
estimate_ar1_noise <- function(y1, y2, trend = NULL) {
  if (length(y1) != length(y2))
    stop("`y1` and `y2` must have the same length.")

  if (is.null(trend)) {
    valid_t <- which(!is.na(y1) & !is.na(y2))
    if (length(valid_t) < 10L)
      stop("Fewer than 10 valid observations; cannot estimate noise.")
    resid1 <- y1[valid_t]
    resid2 <- y2[valid_t]
  } else {
    if (length(y1) != length(trend))
      stop("`trend` must have the same length as `y1` and `y2`.")
    valid_t <- which(!is.na(trend))
    if (length(valid_t) < 10L)
      stop("Fewer than 10 valid trend points; cannot estimate noise.")
    resid1 <- y1[valid_t] - trend[valid_t]
    resid2 <- y2[valid_t] - trend[valid_t]
  }

  list(
    series1 = .variogram_ar1(resid1),
    series2 = .variogram_ar1(resid2)
  )
}
