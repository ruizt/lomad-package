#' Morro Bay dissolved oxygen and pH monitoring data
#'
#' A contiguous block of quality-controlled dissolved oxygen and pH data from
#' the Bay Mouth (BM1) water-quality station in Morro Bay, California,
#' collected through the Central and Northern California Ocean Observing
#' System (CeNCOOS), covering February to April 2022. The test flags one
#' sustained decoupling episode in this block.
#'
#' Raw hourly sensor data were quality-controlled with a QARTOD-based pipeline
#' and split into contiguous blocks; each variable was standardized (z-score);
#' tidal periodicity was removed by zeroing the Fourier components at the five
#' dominant tidal constituents and their spring-neap sidebands; and the series
#' were downsampled to 6-hourly resolution. Values are therefore dimensionless
#' standardized anomalies, not raw concentrations. The processing scripts live
#' in the `lomad-analysis` repository (`mb-analysis/`).
#'
#' @format A data frame with 301 rows and 3 variables:
#' \describe{
#'   \item{datetime}{POSIXct timestamp (6-hourly).}
#'   \item{o2}{Numeric. Standardized, tidally presmoothed dissolved oxygen.}
#'   \item{ph}{Numeric. Standardized, tidally presmoothed pH.}
#' }
#'
#' @source Central and Northern California Ocean Observing System (CeNCOOS)
#'   water-quality monitoring, Morro Bay BM1 station.
#'
#' @seealso [lomad()]
#'
#' @examples
#' plot(morro_bay$datetime, morro_bay$o2, type = "l", col = "steelblue")
#' lines(morro_bay$datetime, morro_bay$ph, col = "tomato")
"morro_bay"
