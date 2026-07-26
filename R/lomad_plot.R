#' Plot detected decoupling periods from a lomad fit
#'
#' Two-panel base R plot. The upper panel shows the MA-smoothed series
#' (`ma1`, `ma2`) and shared trend; the lower panel plots the rolling
#' correlation `R_t` with theoretical `rho_t` as a dashed reference.
#'
#' The two panels shade *different* index ranges, deliberately. A rejection at
#' time \eqn{t} concerns the rolling window
#' \eqn{W_t = \{t - s + 1, \ldots, t\}}, so the lower panel shades the rejected
#' \eqn{t} themselves (where \eqn{R_t} is plotted), while the upper panel
#' shades back to the start of the window, \eqn{t - s + 1}, covering the
#' smoothed values that actually entered the correlation. Upper-panel bands are
#' therefore wider than lower-panel bands by up to \eqn{s - 1} points and are
#' shifted left. The shading stops at the moving averages: the inference is
#' about decoupling of the smoothed series, so the further \eqn{h - 1} raw
#' observations behind each smoothed value are not included.
#'
#' @param fit List returned by [lomad_fit()].
#' @param tst List returned by [lomad_test()].
#' @param dates Optional vector of x-axis values, length `fit$inputs$n`. If it
#'   carries a `Date` or `POSIXct` class the axis is drawn as a calendar axis;
#'   any other numeric vector is used as-is.
#' @param alpha Numeric. Fill transparency for shaded regions (default 0.25).
#' @param band Logical. If `TRUE` (default), shade the null tolerance band in
#'   the lower panel: the region between \eqn{\rho_t} and the critical value
#'   \eqn{\rho_t + z_{\alpha_{\mathrm{eff}}}\sqrt{V_t/s}} below which
#'   \eqn{R_t} is flagged. Because \eqn{\rho_t} and \eqn{V_t} both vary with
#'   \eqn{t}, the band explains why the deepest dip in \eqn{R_t} is not
#'   necessarily the flagged one.
#'
#' @return Invisibly returns `NULL`. Called for its side effect (base R plot).
#'
#' @export
lomad_plot <- function(fit, tst, dates = NULL, alpha = 0.25,
                       band = TRUE) {

  if (is.null(tst))
    stop("lomad_plot() requires a `tst` argument (lomad_test result).")

  n     <- fit$inputs$n
  s     <- fit$inputs$s
  R     <- fit$R
  rho   <- fit$rho
  V     <- fit$V
  ma1   <- fit$ma1
  ma2   <- fit$ma2
  trend <- fit$trend

  rejected <- tst$rejected
  rejected[is.na(rejected)] <- FALSE

  if (!is.null(dates) && length(dates) != n)
    stop(sprintf("`dates` must have length %d (the series length), not %d.",
                 n, length(dates)))

  t_idx <- if (is.null(dates)) seq_len(n) else dates

  rej_shifted <- .rejected_window_span(rejected, s)
  r_bs   <- rle(rej_shifted)
  ends   <- cumsum(r_bs$lengths)
  starts <- ends - r_bs$lengths + 1L
  shade_starts <- starts[r_bs$values]
  shade_ends   <- ends[r_bs$values]

  # Lower panel shading: no back-shift
  r_lo   <- rle(rejected)
  ends_lo   <- cumsum(r_lo$lengths)
  starts_lo <- ends_lo - r_lo$lengths + 1L
  shade_lo_starts <- starts_lo[r_lo$values]
  shade_lo_ends   <- ends_lo[r_lo$values]

  shade_col <- grDevices::rgb(0.7, 0.85, 1, alpha)
  col_trend <- grDevices::rgb(0.4, 0.4, 0.4, 0.8)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::par(mfrow = c(2, 1), oma = c(3, 0, 0, 0))

  # Upper panel: smoothed series + trend
  graphics::par(mar = c(0, 4, 2, 1))
  ylim_top <- range(c(ma1, ma2), na.rm = TRUE)
  ylim_top <- ylim_top + c(-1, 1) * diff(ylim_top) * 0.05
  graphics::plot(t_idx, ma1, type = "n", ylim = ylim_top,
                 xlab = "", ylab = "smoothed series", xaxt = "n")
  .shade_intervals(t_idx, shade_starts, shade_ends, shade_col)
  graphics::lines(t_idx, ma1, col = "blue", lwd = 1.5)
  graphics::lines(t_idx, ma2, col = "red", lwd = 1.5)
  graphics::lines(t_idx, trend, col = col_trend, lwd = 1.2)

  # Lower panel: R_t with rho_t reference
  graphics::par(mar = c(0, 4, 0, 1))
  # Null tolerance band: the test rejects exactly where R_t falls below
  # crit_t = rho_t + z_{alpha_eff} sqrt(V_t / s). Drawing it explains why the
  # deepest dip in R_t need not be significant — rho_t and V_t move too.
  crit <- if (band && !is.null(tst$alpha_eff))
    rho + stats::qnorm(tst$alpha_eff) * sqrt(V / s) else NULL

  R_range  <- range(c(R, rho, crit), na.rm = TRUE)
  ylim_bot <- c(min(R_range[1], -0.1) - 0.05, max(R_range[2], 0.1) + 0.05)
  graphics::plot(t_idx, R, type = "n", ylim = ylim_bot,
                 xlab = "", ylab = "correlation", xaxt = "n")
  .shade_intervals(t_idx, shade_lo_starts, shade_lo_ends, shade_col)
  if (!is.null(crit))
    .band_polygon(t_idx, crit, rho, grDevices::rgb(0.55, 0.55, 0.55, 0.22))
  graphics::lines(t_idx, R, col = "grey40")
  graphics::lines(t_idx, rho, col = "grey30", lty = 2, lwd = 1)
  graphics::abline(h = 0, col = "grey80", lty = 1, lwd = 0.5)

  # Draw a calendar axis when `dates` carries a date/time class; plotting
  # against POSIXct or Date otherwise labels the axis in seconds or days
  # since the epoch.
  if (inherits(t_idx, "POSIXt")) {
    graphics::axis.POSIXct(1, x = t_idx)
  } else if (inherits(t_idx, "Date")) {
    graphics::axis.Date(1, x = t_idx)
  } else {
    graphics::axis(1)
  }

  invisible(NULL)
}


# Shade rectangular regions between start/end index pairs.
.shade_intervals <- function(t_idx, starts, ends, col) {
  if (length(starts) == 0L) return(invisible(NULL))
  usr <- graphics::par("usr")
  for (k in seq_along(starts)) {
    graphics::rect(xleft = t_idx[starts[k]], ybottom = usr[3],
                   xright = t_idx[ends[k]], ytop = usr[4],
                   col = col, border = NA)
  }
}


# Expand a logical vector of rejections into the rolling windows they concern.
#
# A rejection at time t is evidence about W_t = {t - s + 1, ..., t}, so the
# upper panel shades back to the start of that window rather than marking t
# alone. The shift is s - 1 (the correlation window), not h - 1: the panel
# draws the moving averages, and the inferential claim is about decoupling of
# those moving averages rather than of the raw series, so the further h - 1
# raw observations behind each smoothed value are not included.
#
# rejected : logical vector (NA treated as FALSE)
# s        : rolling correlation window length
.rejected_window_span <- function(rejected, s) {
  rejected[is.na(rejected)] <- FALSE
  out <- rep(FALSE, length(rejected))
  s   <- as.integer(s)
  for (t in which(rejected)) out[max(1L, t - (s - 1L)):t] <- TRUE
  out
}


# Fill between two curves over contiguous runs where both are finite.
.band_polygon <- function(x, lo, hi, col) {
  ok <- is.finite(lo) & is.finite(hi)
  if (!any(ok)) return(invisible(NULL))
  r      <- rle(ok)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  for (k in which(r$values)) {
    i <- starts[k]:ends[k]
    graphics::polygon(c(x[i], rev(x[i])), c(lo[i], rev(hi[i])),
                      col = col, border = NA)
  }
  invisible(NULL)
}
