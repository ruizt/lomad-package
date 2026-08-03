#' Fit the LOMAD local correlation model
#'
#' Estimates quantities needed for inference on the local correlation between
#' two time series. Runs the CLT pipeline: MA trend estimation, noise
#' fitting, and rolling computation of \eqn{R_t}, \eqn{\rho_t},
#' and \eqn{V_t}.
#'
#' @param x1 Numeric vector. First time series.
#' @param x2 Numeric vector. Second time series, same length as `x1`.
#' @param h Integer or NULL. Smoothing window half-bandwidth. Auto-selected
#'   if `NULL`.
#' @param s Integer or NULL. Rolling correlation window length.
#'   Auto-selected if `NULL`.
#' @param noise_method Character. Noise estimation method. Only `"ar1"` is
#'   provided: a variogram-based AR(1) fit. Ignored when `noise_override` is
#'   supplied. To use any other noise model, estimate its autocovariance
#'   yourself and pass it through `noise_override`; see `vignette("lomad")`
#'   for a worked example.
#' @param noise_override Optional list giving the noise directly instead of
#'   estimating it from data: either parametric (elements `ar`, `ma`
#'   (optional), `sigma2`) or a raw autocovariance sequence (element `acov`
#'   holding \eqn{\gamma(0), \gamma(1), \ldots}; zero beyond its length).
#'   Useful for oracle experiments or externally estimated noise. May also
#'   be a list of two such lists (one per series).
#' @param min_lambda Non-negative numeric. Windows whose smaller
#'   signal-to-noise ratio falls below this are marked untestable and excluded
#'   from `valid_idx`. Defaults to 0, which excludes nothing and reproduces
#'   the behaviour of earlier versions. Raising it removes windows the test
#'   cannot speak to, but note the resulting subset is not independent of the
#'   p-values: discarded windows are conservative under the null, so filtering
#'   consumes conservatism rather than being free.
#' @param lag_max Integer. Lag truncation for autocovariance sums.
#'   Default 100.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{method}{Character `"clt"`.}
#'     \item{ma1, ma2}{MA-smoothed series.}
#'     \item{noise}{Per-series noise parameter estimates.}
#'     \item{tau1_sq, tau2_sq, rho, V, R}{CLT quantities for inference.}
#'     \item{r_hat}{Affine-invariant effect-size estimate \eqn{R_t/\rho_t^{(0)}}.}
#'     \item{lambda1, lambda2}{Per-window signal-to-noise \eqn{\tau_k^2/\sigma_k^2}.}
#'     \item{testable}{`FALSE` where `min(lambda1, lambda2) < min_lambda`. The
#'       test has no power as either \eqn{\lambda_k \to 0}: \eqn{\rho^{(0)}}
#'       falls to zero with it, so \eqn{R} has nothing to fall below. Such a
#'       window is untestable, which is not the same as showing no decoupling.}
#'     \item{valid_idx}{Integer indices of time points with valid test stats.}
#'     \item{inputs}{List of input parameters (n, h, s, lag_max).}
#'   }
#'
#' @references
#' Ruiz, T. D., Seifert, A. J., Hamilton, E., Mispagel, C. M., Hunt, O. P.,
#' Garcia, J., and Bockmon, E. E. (2026). Inference for local trend similarity
#' in nonstationary time series via rolling correlation, with application to
#' assessing stability in an estuarine system. Manuscript in preparation.
#'
#' @seealso [lomad_test()], [lomad()], [lomad_plot()]
#'
#' @examples
#' # Coupled scenario: the null model holds everywhere, so R_t tracks rho_t
#' tr  <- sim_trends(500, d = 0, method = "smooth", bw = 50, seed = 101)
#' sim <- suppressMessages(
#'   sim_noise_pair(tr, h = 5, lambda_target = 1.5, ar.coefs = 0.5, seed = 102))
#' fit <- lomad_fit(sim$y1, sim$y2, h = 5, s = 125)
#' plot(fit$R, type = "l", ylab = "R_t")
#' lines(fit$rho, lwd = 2)               # null benchmark rho_t
#'
#' # Real data at the paper's settings
#' fit_mb <- lomad_fit(morro_bay$o2, morro_bay$ph, h = 4, s = 60)
#'
#' @export
lomad_fit <- function(x1             = NULL,
                      x2             = NULL,
                      h              = NULL,
                      s              = NULL,
                      noise_method   = c("ar1"),
                      noise_override = NULL,
                      lag_max        = 100L,
                      min_lambda     = 0) {

  if (is.null(x1) || is.null(x2))
    stop("`x1` and `x2` are required.")
  if (length(x1) != length(x2))
    stop("`x1` and `x2` must have the same length.")

  noise_method <- match.arg(noise_method)

  .lomad_fit_clt(y1 = x1, y2 = x2, h = h, s = s, lag_max = lag_max,
                 min_lambda = min_lambda,
                 noise_method = noise_method, noise_override = noise_override)
}
