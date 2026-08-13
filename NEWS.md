# lomad 0.1.0

Implements the two-tau form of Proposition 1: each series carries its own
signal variance, and the shared trend is gone from the API.

## Breaking changes

* `estimate_ar1_noise()` takes `trend1`, `trend2` and `h` in place of a single
  `trend`. **A 0.0.1 positional call `estimate_ar1_noise(y1, y2, trend)` still
  runs, but detrends series 1 only** and leaves series 2 raw. Pass the same
  vector as both `trend1` and `trend2` for the previous behaviour.

* `compute_rho()` and `compute_V()` take `tau1_sq` and `tau2_sq` in place of a
  single `tau_sq`.

* `lomad_fit()` no longer returns `trend`. Use `ma1` and `ma2`.

* `lomad_plot()` no longer draws the shared trend.

## New

* `lomad_fit()` gains `min_lambda`, defaulting to `0`, and returns `testable`,
  `r_hat`, `lambda1` and `lambda2`.

* `sim_trends()` gains an optional affine layer, `x2 <- a_t + b_t * x2`, via
  `affine_s`, `affine_bw` and `affine_cap`. The coefficients drift slowly
  enough to be near constant within any window of length `affine_s` while
  accumulating across the series, so the pair stays locally affine similar but
  distant windows see different maps. It is applied after `d`, which continues
  to act on the base trends as before, and is off by default. The return value
  gains the realized `a` and `b`.

  The coefficients can instead be supplied outright via `affine_a` and
  `affine_b`, and `affine_cap = 0` gives a drift-free layer. Neither draws any
  random numbers, so both leave a seed behaving exactly as it did before this
  argument existed: `affine_a = rep(0, n)` with `affine_b = rep(2, n)`
  reproduces a fixed `x2 <- 2 * x2` bit for bit.

## Changed

* `estimate_ar1_noise()` estimates noise from per-series residuals rather than
  from residuals against a shared trend, and corrects the AR(1) estimate for
  the high-pass filtering that subtracting a moving average induces when `h` is
  supplied.

* `sim_trends()` normalises the first coefficient vector to fixed total power
  and constrains the displacement to be orthogonal to it, rescaled to equal
  norm. Both act on the base coefficient pair shared by every `method`, so a
  given `seed` produces different trends than it did in 0.0.1.
