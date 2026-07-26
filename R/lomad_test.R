#' Test for local correlation decoupling
#'
#' Given the output of [lomad_fit()], tests whether the local sample
#' correlation is significantly below the estimated population correlation
#' using the BY-corrected pointwise Z-test.
#'
#' @param fit List returned by [lomad_fit()].
#' @param alpha Numeric. FDR level for BY correction (default 0.05).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{Z}{Numeric vector of Z-statistics (NA outside valid indices).}
#'     \item{p_values}{Raw p-values.}
#'     \item{p_adj}{BY-adjusted p-values.}
#'     \item{rejected}{Logical vector of rejections.}
#'     \item{alpha_eff}{Effective significance level after BY correction.}
#'   }
#'
#' @references
#' Ruiz, T. D., Seifert, A. J., Hamilton, E., Mispagel, C. M., Hunt, O. P.,
#' Garcia, J., and Bockmon, E. E. (2026). Inference for local trend similarity
#' in nonstationary time series via rolling correlation, with application to
#' assessing stability in an estuarine system. Manuscript in preparation.
#'
#' Benjamini, Y. and Yekutieli, D. (2001). The control of the false discovery
#' rate in multiple testing under dependency. \emph{The Annals of Statistics},
#' 29(4), 1165--1188.
#'
#' @seealso [lomad_fit()], [lomad()]
#'
#' @examples
#' dcp <- subset(sim_decoupling, scenario == "decoupled")
#' fit <- lomad_fit(dcp$y1, dcp$y2, h = 5, s = 125)
#' tst <- lomad_test(fit, alpha = 0.05)
#' range(which(tst$rejected))            # extent of flagged region
#' tst$alpha_eff                         # BY-adjusted threshold
#'
#' @export
lomad_test <- function(fit, alpha = 0.05) {
  .lomad_test_clt(fit, alpha = alpha)
}
