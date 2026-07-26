# data-raw/morro_bay.R
#
# Packages the `morro_bay` example dataset: two contiguous blocks of
# quality-controlled dissolved oxygen and pH monitoring data from the Bay
# South (BS1) station in Morro Bay, California (CeNCOOS water-quality
# sensors).
#
# Provenance (see the lomad-analysis repository, mb-analysis/):
#   1. Raw 10-15 min sensor data quality-controlled with a QARTOD-based
#      pipeline (mb-qartod repository) and split into contiguous blocks.
#   2. Each variable standardized (z-score).
#   3. Tidal periodicity removed with a centered 25 h rolling mean, then
#      downsampled to 6-hourly resolution.
#
# Source file: lomad-analysis/example-data/mb-example.rds (sibling repo).
#
# Run from the package root: source("data-raw/morro_bay.R")

src <- "../lomad-analysis/example-data/mb-example.rds"
if (!file.exists(src))
  stop("mb-example.rds not found; expected lomad-analysis as a sibling repo.")

mb <- readRDS(src)

morro_bay <- do.call(rbind, lapply(mb, function(df) {
  data.frame(
    block    = as.integer(df$block_id),
    datetime = df$datetime,
    o2       = df$o2,
    ph       = df$ph
  )
}))
rownames(morro_bay) <- NULL

save(morro_bay, file = "data/morro_bay.rda", compress = "xz")
cat("Wrote data/morro_bay.rda:", nrow(morro_bay), "rows,",
    length(unique(morro_bay$block)), "blocks\n")
