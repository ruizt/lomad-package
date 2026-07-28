# data-raw/morro_bay.R
#
# Packages the `morro_bay` example dataset from the static subset exported by
# the analysis repository.
#
# Input:  data-raw/mb-example.rds   (committed; no cross-repo dependency)
# Output: data/morro_bay.rda
#
# Run from the package root:
#   source("data-raw/morro_bay.R")
#
# The .rds is produced by lomad-analysis/mb-analysis/export_example.R and
# handed over as a file. That is deliberately one-way: this package never reads
# from the analysis repository, so `morro_bay` can be rebuilt from a clean
# checkout. Re-export only when the example blocks should change.
#
# Provenance of the subset (all upstream of the hand-off):
#   1. Raw hourly sensor data quality-controlled with a QARTOD-based pipeline
#      (mb-qartod).
#   2. Binned hourly, split into contiguous blocks at gaps over 24 h, blocks
#      under 5 days dropped, each variable standardized globally
#      (lomad-analysis, mb-analysis/process_blocks.R).
#   3. Tidal periodicity removed by zeroing the Fourier components at the five
#      dominant tidal constituents and their spring-neap sidebands, then
#      downsampled to 6-hourly (lomad-analysis, mb-analysis/utils.R).
#
# Step 3 previously used a centered 25 h rolling mean. That is a low-pass whose
# attenuation reaches well below its nominal cutoff, so the 6-hourly series it
# produced was band-limited far below its own Nyquist frequency: the lag-2
# variogram exceeded twice the lag-1 variogram, which no stationary AR(1) can
# produce, and lomad_fit() clamped phi_hat at its boundary on every fit.

src <- "data-raw/mb-example.rds"
if (!file.exists(src))
  stop(src, " not found. It is exported by lomad-analysis ",
       "(mb-analysis/export_example.R) and committed to this repository.")

morro_bay <- readRDS(src)

stopifnot(is.data.frame(morro_bay),
          identical(names(morro_bay), c("block", "datetime", "o2", "ph")),
          nrow(morro_bay) > 0L,
          !anyNA(morro_bay))

save(morro_bay, file = "data/morro_bay.rda", compress = "xz")
cat(sprintf("Wrote data/morro_bay.rda: %d rows, %d blocks (%s)\n",
            nrow(morro_bay), length(unique(morro_bay$block)),
            paste(unique(morro_bay$block), collapse = ", ")))

# The shipped data must be something the noise model can actually fit.
if (requireNamespace("lomad", quietly = TRUE)) {
  vg <- function(x, l) mean((x[(l + 1):length(x)] - x[1:(length(x) - l)])^2) / 2
  for (b in unique(morro_bay$block)) {
    d <- morro_bay[morro_bay$block == b, ]
    f <- suppressWarnings(suppressMessages(
      lomad::lomad_fit(d$o2, d$ph, h = 4, s = 60)))
    vt <- which(!is.na(f$trend)); r <- d$o2[vt] - f$trend[vt]
    cat(sprintf("  block %2d: n = %3d, V2/V1 = %.2f, phi_hat = %.3f\n",
                b, nrow(d), vg(r, 2) / vg(r, 1), f$noise$series1$ar))
  }
  cat("  (V2/V1 must be below 2 for a stationary AR(1) to fit)\n")
}
