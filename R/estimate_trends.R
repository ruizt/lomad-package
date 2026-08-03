#' Estimate per-series trends by MA smoothing
#'
#' Applies a one-sided moving-average filter of width \eqn{h} to each series:
#' \deqn{\hat\nu^{(h)}_{kt} = \frac{1}{h}\sum_{j=0}^{h-1} y_{k,t-j}.}
#' Under affine similarity the two trends need not coincide, so no shared
#' trend is formed.
#' The first \eqn{h - 1} observations are `NA` due to the one-sided filter.
#'
#' @param y1 Numeric vector. First observed series.
#' @param y2 Numeric vector. Second observed series, same length as `y1`.
#' @param h Positive integer. MA smoothing window width (paper notation:
#'   \eqn{h}).
#'
#' @return A named list:
#'   \describe{
#'     \item{ma1}{Numeric vector (length `n`). MA-smoothed series 1.}
#'     \item{ma2}{Numeric vector (length `n`). MA-smoothed series 2.}
#'     \item{h}{Integer. Window width used.}
#'   }
#'
#' @export
estimate_trends <- function(y1, y2, h) {
  if (length(y1) != length(y2))
    stop("`y1` and `y2` must have the same length.")
  h  <- as.integer(h)
  if (h < 1L) stop("`h` must be >= 1.")

  kernel <- rep(1 / h, h)
  ma1    <- as.numeric(stats::filter(y1, kernel, sides = 1))
  ma2    <- as.numeric(stats::filter(y2, kernel, sides = 1))

  list(ma1 = ma1, ma2 = ma2, h = h)
}
