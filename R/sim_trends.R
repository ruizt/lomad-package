#' Simulate a pair of Fourier-basis trend series with controlled separation
#'
#' Generates two time series from a shared Fourier basis. The coupling weight
#' `w` controls the time-varying mixing between the two underlying trends:
#' `w = 0` means fully decoupled (maximum separation), `w = 1` means fully
#' coupled (both series track their shared mean). `d` scales the distinct
#' component against that shared mean, so separation grows linearly in `d`
#' while the coupling profile is left as the structure produced it.
#'
#' An optional affine layer may then be applied to the second series,
#' `x2 <- a_t + b_t * x2`, with `a_t` and `b_t` drifting slowly enough that they
#' are near constant within any window of length `affine_s`. Local affine
#' similarity constrains the coefficients only within a window, so capping the
#' per-window increment leaves them free to accumulate across the series: the
#' null still holds everywhere, while windows far apart see genuinely different
#' maps. The layer is applied after `d` has been imposed, so `d` continues to
#' operate on the base trends exactly as it does without it. The coefficients
#' can also be supplied outright via `affine_a` and `affine_b`, which draws no
#' random numbers and so leaves any seed behaving as it did before the layer
#' existed.
#'
#' The coupling weight can be specified in three ways:
#' \enumerate{
#'   \item **Named method** (`method`): one of `"dist"` (static separation,
#'     default), `"rs"` (stochastic repulsion), `"rm"` (stochastic crossing),
#'     or `"fr"` (fixed-rate, event-based decoupling). Method-specific
#'     parameters are passed via `...`.
#'   \item **Numeric vector** (`w`): a precomputed coupling weight of length
#'     `n`, overriding `method`.
#'   \item **Function** (`w`): a function `f(n)` returning a numeric vector of
#'     length `n`, called at generation time.
#' }
#'
#' @param n Integer. Length of the output series (default 500).
#' @param d Numeric. Amplitude of the distinct component relative to the shared
#'   mean trend (default 1). Separation is linear in `d`, but `d` is not the
#'   realized separation: how much local separation a given `d` produces
#'   depends on the coupling structure and on the window it is measured over.
#'   The realized quantity is what the study reports.
#' @param method Character. Coupling method when `w` is not supplied. One of
#'   `"dist"` (default), `"rs"`, `"rm"`, or `"fr"`.
#' @param w Coupling weight: a numeric vector of length `n`, a function
#'   `f(n) -> numeric(n)`, or `NULL` (default, uses `method`).
#' @param nb Integer. Number of Fourier basis functions (must be odd,
#'   default 25).
#' @param sd0 Numeric. Standard deviation of the lowest-frequency Fourier
#'   coefficient (default 2). Higher-frequency coefficients decay as
#'   `sd0 / k^p`.
#' @param p Numeric. Spectral decay exponent (default 2.5).
#' @param k_min Integer. Minimum harmonic index to include (default 1). Setting
#'   `k_min > 1` excludes low-frequency components, concentrating signal power
#'   at shorter periods.
#' @param affine_a,affine_b Affine coefficients supplied directly: a numeric
#'   vector of length `n`, a function `f(n) -> numeric(n)`, or `NULL`
#'   (default). Supplying either one overrides `affine_s` and takes the
#'   coefficients literally, on the scale of the data. Whichever is `NULL`
#'   falls back to its identity value (`0` for `affine_a`, `1` for
#'   `affine_b`). Nothing is drawn in this case, so a fixed map applied to a
#'   seed predating these arguments is recovered exactly: `affine_a = rep(0, n)`
#'   with `affine_b = rep(2, n)` reproduces `x2 <- 2 * x2`.
#' @param affine_s Integer or NULL. Window length the affine drift cap is
#'   calibrated to, when the coefficients are generated rather than supplied.
#'   `NULL` (default) disables the affine layer, leaving `a_t = 0` and
#'   `b_t = 1`. Set it to the window length the series will be analysed at.
#' @param affine_bw Numeric. Smoothing bandwidth of the coefficient paths, as a
#'   fraction of `n` (default 0.5). Larger values give smoother paths that
#'   accumulate more total drift from the same per-window budget; smaller ones
#'   turn more often and accumulate less.
#' @param affine_cap Numeric. Maximum drift in `a_t` and `b_t` over any window
#'   of length `affine_s` (default 0.015), as a fraction of a series standard
#'   deviation for `a_t` and of unit slope for `b_t`. `0` gives a drift-free
#'   layer, identical in output and in RNG consumption to leaving `affine_s`
#'   unset, so a sweep can carry the no-drift arm as one more cap value.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#' @param ... Additional arguments passed to the coupling weight generator when
#'   using a named `method`:
#'   \describe{
#'     \item{`"rs"`}{`bw` (bandwidth, default 50), `coupling` (fraction of
#'       time in coupled state, default 0.8).}
#'     \item{`"rm"`}{`bw` (bandwidth, default 50), `coupling` (fraction of
#'       time near coupled state, default 0.8).}
#'     \item{`"fr"`}{`rate` (decoupling events per unit time, default 0.01)
#'       and `bump`, the pulse shape: `"gaussian"` (default) or `"gamma"`.
#'       Both place identical events at identical times and differ only in
#'       smoothness at onset -- the gamma pulse has a corner there, the
#'       gaussian does not. Because difference-based noise estimation is
#'       trend-robust only for trends with bounded derivative (Hall and
#'       Van Keilegom, 2003), the two behave very differently under strong
#'       noise autocorrelation, which is why the smooth shape is the default.}
#'   }
#'
#' @return A list with:
#'   \describe{
#'     \item{x1}{Numeric vector of length `n`. First output series.}
#'     \item{x2}{Numeric vector of length `n`. Second output series, after the
#'       affine layer if one was applied.}
#'     \item{x_mean}{Numeric vector of length `n`. Shared mean trend of the base
#'       trends, before any affine layer.}
#'     \item{w}{Numeric vector of length `n`. Coupling weight used.}
#'     \item{a}{Numeric vector of length `n`. Realized affine intercept, all
#'       zero when the layer is disabled.}
#'     \item{b}{Numeric vector of length `n`. Realized affine gradient, all one
#'       when the layer is disabled.}
#'   }
#'
#' @examples
#' # Unstructured (static separation)
#' tr <- sim_trends(500, d = 2)
#'
#' # Named method with tuning arguments
#' tr <- sim_trends(500, d = 2, method = "rs", bw = 50, coupling = 0.8)
#'
#' # Custom coupling weight vector
#' w_custom <- rep(c(0, 1), each = 250)
#' tr <- sim_trends(500, d = 2, w = w_custom)
#'
#' # Locally affine-similar pair: identical up to a slowly drifting map
#' tr <- sim_trends(2500, d = 0, method = "rs", affine_s = 100)
#' range(tr$b)
#'
#' # A fixed map, supplied rather than generated
#' tr <- sim_trends(500, d = 0, affine_b = rep(2, 500))
#'
#' @references
#' Hall, P. and Van Keilegom, I. (2003). Using difference-based methods for
#' inference in nonparametric regression with time series errors. \emph{Journal
#' of the Royal Statistical Society Series B}, 65(2), 443--456.
#'
#' @export
sim_trends <- function(n          = 500,
                       d          = 1,
                       method     = c("dist", "rs", "rm", "fr"),
                       w          = NULL,
                       nb         = 25,
                       sd0        = 2,
                       p          = 2.5,
                       k_min      = 1L,
                       affine_a   = NULL,
                       affine_b   = NULL,
                       affine_s   = NULL,
                       affine_bw  = 0.5,
                       affine_cap = 0.015,
                       seed       = NULL,
                       ...) {
  if ((nb %% 2) == 0) stop("`nb` must be odd.")

  given <- !is.null(affine_a) || !is.null(affine_b)
  if (!is.null(affine_s)) {
    affine_s <- as.integer(affine_s)
    if (is.na(affine_s) || affine_s < 2L || affine_s > n)
      stop("`affine_s` must be between 2 and `n` (", n, ").")
    if (affine_bw <= 0) stop("`affine_bw` must be positive.")
    if (affine_cap < 0) stop("`affine_cap` must be non-negative.")
  }
  resolve_coef <- function(x, default, nm) {
    if (is.null(x)) return(rep(default, n))
    if (is.function(x)) x <- x(n)
    if (!is.numeric(x) || length(x) != n)
      stop("`", nm, "` must be numeric of length `n` (", n, "), got ",
           length(x), ".")
    x
  }

  # Resolve coupling weight
  if (is.null(w)) {
    method <- match.arg(method)
    w <- switch(method,
      dist = rep(0, n),
      rs   = .make_w_smooth(n, ...),
      rm   = .make_w_cross(n, ...),
      fr   = .make_w_rate(n, ...)
    )
  } else if (is.function(w)) {
    w <- w(n)
  }

  if (length(w) != n) {
    stop("`w` must have length `n` (", n, "), got ", length(w), ".")
  }

  # Build base trends from Fourier basis
  basis <- .make_basis_trends(n = n, nb = nb, sd0 = sd0, d = 1,
                              p = p, k_min = k_min, seed = seed)

  # Mix and rescale
  out <- .apply_w(basis$mu1, basis$mu2, basis$x_mean, w, d)

  # Affine layer last, so `d` is imposed on the base trends either way.
  #
  # Supplied coefficients draw nothing, which is what lets a fixed map be
  # recovered from a seed that predates this argument.
  if (given) {
    out$a <- resolve_coef(affine_a, 0, "affine_a")
    out$b <- resolve_coef(affine_b, 1, "affine_b")
    out$x2 <- out$a + out$b * out$x2
  } else if (is.null(affine_s)) {
    out$a <- rep(0, n)
    out$b <- rep(1, n)
  } else {
    af     <- .apply_affine(out$x1, out$x2, s = affine_s,
                            bw = affine_bw, cap = affine_cap)
    out$x2 <- af$x2
    out$a  <- af$a
    out$b  <- af$b
  }
  out
}
