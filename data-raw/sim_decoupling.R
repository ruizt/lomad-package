# data-raw/sim_decoupling.R
#
# Generates the packaged `sim_decoupling` example dataset: two simulated
# scenarios sharing the same base trends and noise process, differing only in
# the L2 trend separation d.
#
#   coupled:    d = 0    (exactly shared trend; the null hypothesis holds)
#   decoupled:  d = 1.5  (structured separation via smooth coupling weight)
#
# DGP mirrors the paper's power study: sim_trends(method = "smooth", bw = 50)
# + sim_noise_pair() with AR(1) noise (phi = 0.5) calibrated to SNR 1.5 at
# smoothing bandwidth h = 5.
#
# Run from the package root: source("data-raw/sim_decoupling.R")

devtools::load_all(".")

n    <- 500L
h    <- 5L     # smoothing bandwidth used for SNR calibration (and analysis)
phi  <- 0.5    # AR(1) noise coefficient
snr  <- 1.5    # target signal-to-noise ratio
seed <- 101L

make_scenario <- function(d, label) {
  trends <- sim_trends(n = n, d = d, method = "smooth", bw = 50, seed = seed)
  sim    <- sim_noise_pair(trends, h = h, lambda_target = snr,
                           ar.coefs = phi, seed = seed + 1L)
  data.frame(
    scenario = label,
    t        = seq_len(n),
    y1       = sim$y1,
    y2       = sim$y2,
    nu1      = trends$x1,
    nu2      = trends$x2
  )
}

sim_decoupling <- rbind(
  make_scenario(d = 0,   label = "coupled"),
  make_scenario(d = 1.5, label = "decoupled")
)
sim_decoupling$scenario <- factor(sim_decoupling$scenario,
                                  levels = c("coupled", "decoupled"))

save(sim_decoupling, file = "data/sim_decoupling.rda", compress = "xz")
cat("Wrote data/sim_decoupling.rda:", nrow(sim_decoupling), "rows\n")
