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

## Changed

* `estimate_ar1_noise()` estimates noise from per-series residuals rather than
  from residuals against a shared trend, and corrects the AR(1) estimate for
  the high-pass filtering that subtracting a moving average induces when `h` is
  supplied.

* `sim_trends()` normalises the first coefficient vector to fixed total power
  and constrains the displacement to be orthogonal to it, rescaled to equal
  norm. Both act on the base coefficient pair shared by every `method`, so a
  given `seed` produces different trends than it did in 0.0.1.
