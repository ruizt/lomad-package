# data-raw/morro_bay.R
#
# Packages the `morro_bay` example dataset: two contiguous blocks of
# quality-controlled dissolved oxygen and pH monitoring data from the Bay
# South (BS1) station in Morro Bay, California (CeNCOOS water-quality
# sensors).
#
# Provenance — the full chain, reproducible from sibling repositories:
#   1. Raw hourly sensor data quality-controlled with a QARTOD-based pipeline
#      (mb-qartod, inst/scripts/run_qartod_pipeline.R -> wp_data.parquet).
#   2. Binned hourly, split into contiguous blocks at gaps over 24 h, blocks
#      under 5 days dropped, each variable standardized globally
#      (lomad-analysis, mb-analysis/process_blocks.R -> ph_o2_blocks.csv).
#   3. Tidal periodicity removed by zeroing the Fourier components at the five
#      dominant tidal constituents and their spring-neap sidebands, then
#      downsampled to 6-hourly (lomad-analysis, mb-analysis/utils.R ::
#      presmooth_tidal).
#
# Step 3 previously used a centered 25 h rolling mean. That is a low-pass
# filter, and a boxcar's attenuation reaches well below its nominal cutoff, so
# the 6-hourly series it produced was band-limited far below its own Nyquist
# frequency. The consequence was visible in the shipped data: the lag-2
# variogram exceeded twice the lag-1 variogram, which no stationary AR(1)
# process can produce, so lomad_fit() could not estimate the noise model and
# clamped phi_hat at its boundary on every fit. The notch removes only the
# tidal bands and leaves the rest of the spectrum intact.
#
# Run from the package root, with lomad-analysis as a sibling:
#   source("data-raw/morro_bay.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

analysis   <- "../lomad-analysis"
blocks_csv <- file.path(analysis, "_mb-data/ph_o2_blocks.csv")
utils_r    <- file.path(analysis, "mb-analysis/utils.R")

for (p in c(blocks_csv, utils_r))
  if (!file.exists(p))
    stop(basename(p), " not found at ", p,
         "\nExpected lomad-analysis as a sibling repo with ",
         "mb-analysis/process_blocks.R already run.")

source(utils_r)   # presmooth_tidal(), TIDAL_PERIODS

# Blocks are identified by the period they cover, not by block_id: the ids are
# assigned by position within a station and shift whenever the record is
# extended or the QA changes.
#
# Two consecutive Bay Mouth blocks are shipped: winter 2021-22, where the test
# flags nothing, and the spring 2022 block immediately after it, which contains
# one sustained decoupling episode. The station changed from Bay South with the
# tidal-filter correction. Under the corrected pipeline Bay South records no
# detections in any of its seven analysable blocks, so a Bay South example
# would demonstrate the method finding nothing; Bay Mouth flags episodes in
# seven of nine.
KEEP_LOCATION <- "BM1"
PERIODS <- list(c("2021-10-01", "2022-01-31"),
                c("2022-02-01", "2022-04-30"))
MIN_HOURS <- 1000L   # excludes short fragments that share a window

all_blocks <- read_csv(blocks_csv, show_col_types = FALSE) |>
  filter(location == KEEP_LOCATION)

pick <- vapply(PERIODS, function(w) {
  hit <- all_blocks |>
    group_by(block_id) |>
    summarise(n = n(), start = min(datetime), end = max(datetime),
              .groups = "drop") |>
    filter(start >= as.POSIXct(w[1], tz = "UTC"),
           end   <= as.POSIXct(w[2], tz = "UTC"),
           n >= MIN_HOURS)
  if (nrow(hit) != 1L)
    stop("expected exactly one ", KEEP_LOCATION, " block within ",
         w[1], " .. ", w[2], "; found ", nrow(hit))
  hit$block_id
}, numeric(1))

raw <- all_blocks |> filter(block_id %in% pick)

morro_bay <- raw |>
  group_split(block_id) |>
  lapply(presmooth_tidal, cols = c("o2", "ph"), step = "6 hours") |>
  bind_rows() |>
  transmute(block    = as.integer(block_id),
            datetime = datetime,
            o2       = o2,
            ph       = ph) |>
  arrange(block, datetime) |>
  as.data.frame()
rownames(morro_bay) <- NULL

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
    cat(sprintf("  block %2d: n = %4d, V2/V1 = %.2f, phi_hat = %.3f\n",
                b, nrow(d), vg(r, 2) / vg(r, 1), f$noise$series1$ar))
  }
  cat("  (V2/V1 must be below 2 for a stationary AR(1) to fit)\n")
}
