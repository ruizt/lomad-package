# Internal test implementation for lomad_test()
#
# .lomad_test_clt() — pointwise CLT test (paper method)


# Pointwise CLT test: the paper method.
.lomad_test_clt <- function(fit, alpha) {
  if (alpha <= 0 || alpha >= 1)
    stop("`alpha` must be in (0, 1).")

  n         <- fit$inputs$n
  s         <- fit$inputs$s
  valid_idx <- fit$valid_idx
  m         <- length(valid_idx)

  if (m == 0L) stop("No valid time points in fit.")

  R   <- fit$R
  rho <- fit$rho
  V   <- fit$V

  Z      <- rep(NA_real_, n)
  p_raw  <- rep(NA_real_, n)

  se                  <- sqrt(V[valid_idx] / s)
  Z[valid_idx]        <- (R[valid_idx] - rho[valid_idx]) / se
  p_raw[valid_idx]    <- stats::pnorm(Z[valid_idx])

  # Benjamini-Yekutieli step-up. An earlier version applied the single fixed
  # threshold alpha / c(m) -- the BY boundary evaluated at its largest point,
  # k = m -- to every test. That rejects at least as much as the step-up
  # procedure and controls nothing: under the global null with independent
  # p-values its familywise error grows with m (0.16 at m = 10, ~1 at
  # m = 1000). Heavy window overlap masked this at the m/s ~ 3 used in the
  # simulation studies. The step-up procedure controls FDR at alpha under
  # arbitrary dependence (Benjamini and Yekutieli, 2001), so its validity
  # needs no per-configuration calibration.
  p_adj            <- rep(NA_real_, n)
  p_adj[valid_idx] <- stats::p.adjust(p_raw[valid_idx], method = "BY")

  rejected            <- rep(NA, n)
  rejected[valid_idx] <- p_adj[valid_idx] <= alpha

  # Realized rejection threshold on the raw p-value scale. For a step-up
  # procedure the rejection set is exactly {p_raw <= alpha * k / (m c(m))}
  # with k the number of rejections, so this threshold reproduces `rejected`
  # and gives plots a line consistent with the flags. With k = 0 it falls
  # back to the k = 1 boundary: the bar any single p-value needed to clear.
  c_m        <- sum(1 / seq_len(m))
  n_rejected <- sum(rejected, na.rm = TRUE)
  alpha_eff  <- alpha * max(n_rejected, 1L) / (m * c_m)

  message(sprintf("Rejected %d / %d time points at BY-FDR = %.2f (alpha_eff = %.4g)",
                  n_rejected, m, alpha, alpha_eff))

  list(
    Z         = Z,
    p_values  = p_raw,
    p_adj     = p_adj,
    rejected  = rejected,
    alpha_eff = alpha_eff,
    inputs    = list(alpha = alpha, s = s)
  )
}
