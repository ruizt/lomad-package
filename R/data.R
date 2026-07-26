#' Simulated coupled and decoupled series pairs
#'
#' Two simulated bivariate time series scenarios sharing the same base trends
#' and noise process, differing only in the L2 trend separation \eqn{d}:
#' `"coupled"` (\eqn{d = 0}; the trends are identical, so the local-similarity
#' null hypothesis holds everywhere) and `"decoupled"` (\eqn{d = 1.5};
#' structured separation concentrated in local episodes via a smooth coupling
#' weight).
#'
#' Generated with [sim_trends()] (`method = "smooth"`, `bw = 50`, seed 101)
#' and [sim_noise_pair()] with AR(1) noise (\eqn{\phi = 0.5}) calibrated to a
#' signal-to-noise ratio of 1.5 at smoothing bandwidth \eqn{h = 5}. The
#' generating script is `data-raw/sim_decoupling.R` in the package source
#' repository.
#'
#' @format A data frame with 1000 rows (two scenarios of 500 time points) and
#'   6 variables:
#' \describe{
#'   \item{scenario}{Factor: `"coupled"` or `"decoupled"`.}
#'   \item{t}{Integer time index, 1 to 500.}
#'   \item{y1, y2}{Numeric. Observed series (trend plus AR(1) noise).}
#'   \item{nu1, nu2}{Numeric. True underlying trends (identical in the
#'     coupled scenario).}
#' }
#'
#' @seealso [lomad()], [sim_trends()], [sim_noise_pair()]
#'
#' @examples
#' cpl <- subset(sim_decoupling, scenario == "coupled")
#' plot(cpl$t, cpl$y1, type = "l", col = "grey60")
#' lines(cpl$t, cpl$nu1, lwd = 2)
"sim_decoupling"


#' Morro Bay dissolved oxygen and pH monitoring data
#'
#' Two contiguous blocks of quality-controlled dissolved oxygen and pH data
#' from the Bay South (BS1) water-quality station in Morro Bay, California,
#' collected through the Central and Northern California Ocean Observing
#' System (CeNCOOS). Used as the application example in the accompanying
#' paper.
#'
#' Raw sensor data (10-15 minute resolution) were quality-controlled with a
#' QARTOD-based pipeline and split into contiguous blocks; each variable was
#' standardized (z-score); tidal periodicity was removed with a centered 25 h
#' rolling mean; and the series were downsampled to 6-hourly resolution.
#' Values are therefore dimensionless standardized anomalies, not raw
#' concentrations. The processing scripts live in the `lomad-analysis`
#' repository (`mb-analysis/`).
#'
#' @format A data frame with 1391 rows and 4 variables:
#' \describe{
#'   \item{block}{Integer block identifier (2: Apr-Dec 2020, 971
#'     observations; 21: Sep-Dec 2023, 420 observations).}
#'   \item{datetime}{POSIXct timestamp (6-hourly).}
#'   \item{o2}{Numeric. Standardized, tidally presmoothed dissolved oxygen.}
#'   \item{ph}{Numeric. Standardized, tidally presmoothed pH.}
#' }
#'
#' @source Central and Northern California Ocean Observing System (CeNCOOS)
#'   water-quality monitoring, Morro Bay BS1 station.
#'
#' @seealso [lomad()]
#'
#' @examples
#' b2 <- subset(morro_bay, block == 2)
#' plot(b2$datetime, b2$o2, type = "l", col = "steelblue")
#' lines(b2$datetime, b2$ph, col = "tomato")
"morro_bay"
