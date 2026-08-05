# lomad 0.1.0

The null hypothesis is now trend identity under an arbitrary unknown affine
transformation: \eqn{\nu_2 = a + b\nu_1} with \eqn{b > 0}. Proposition 1 carries
two signal-variance terms rather than one, and the shared-trend approximation is
gone. Results from 0.0.1 will not reproduce under this version.

## Breaking changes

* `estimate_ar1_noise()` takes `trend1`, `trend2` and `h` in place of a single
  `trend`. **A 0.0.1 positional call `estimate_ar1_noise(y1, y2, trend)` still
  runs, but now detrends series 1 only** and leaves series 2 raw. Pass the same
  vector as both `trend1` and `trend2` to reproduce the old behaviour.

* `lomad_fit()` no longer returns `trend`. Under affine similarity the half-sum
  `(ma1 + ma2)/2` estimates nothing when \eqn{b \neq 1}, so the shared trend was
  removed rather than left as a quantity that happens to be defined. Use `ma1`
  and `ma2`.

* `lomad_plot()` no longer draws the shared trend.

* `sim_trends()` returns different trends for a given `seed` than 0.0.1 did, for
  every `method`. Two changes are responsible, both in the base Fourier pair
  that all methods build on: the first coefficient vector is normalised to fixed
  total power, and the displacement is projected orthogonal to it and rescaled
  to equal norm. Seeds are not comparable across versions.

## New

* `lomad_fit()` gains `min_lambda` and returns `testable`, marking windows whose
  smaller local SNR falls below it. Default `0` flags nothing, preserving 0.0.1
  behaviour. Windows with `lambda` near zero have population correlation near
  zero and the test has no power there; `testable` makes that visible rather
  than leaving it to be inferred.

* `lomad_fit()` returns `r_hat`, the affine-invariant effect size
  \eqn{R_t / \rho_t^{(0)}}, and the per-series `lambda1`, `lambda2`.

* `compute_rho()` and `compute_V()` take `tau1_sq` and `tau2_sq` separately.
  Note the cross-pairing in `compute_V()`: \eqn{\tau_2^2} multiplies \eqn{L_1}
  and \eqn{\tau_1^2} multiplies \eqn{L_2}.

## Fixes

* Noise is estimated from per-series residuals. Subtracting a shared trend made
  the \eqn{\hat\phi} bias grow with the separation between the two trends, which
  cost detection where separation was largest.

* `\eqn{\hat\phi}` is corrected for the high-pass filtering that subtracting a
  moving average induces. The residual is an FIR filter of the noise with
  weights \eqn{(1 - 1/h, -1/h, \ldots)}; the variogram ratio is inverted for
  \eqn{\phi} against that filter rather than against the unfiltered process.

## Known limitations

* \eqn{\hat\rho^{(0)}} is biased low by roughly 0.08 wherever there is signal to
  estimate, which makes the test conservative. The bias is the residue of two
  larger opposing errors: `pmax(0, .)` on \eqn{\hat\tau^2} rectifies noise
  upward, and the square root in \eqn{\rho = \sqrt{g(\lambda_1)g(\lambda_2)}}
  pulls it back down. Correcting either alone makes calibration worse.

* Size is inflated near the unit root — about 0.115 against a nominal 0.05 at
  \eqn{\phi = 0.8} — driven by the variance of a single per-series
  \eqn{\hat\phi}. Series with weaker autocorrelation are unaffected.
