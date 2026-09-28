#' Simulated data for illustrating MoCo
#'
#' A simulated data set with 400 observations and seven continuous outcomes.
#' It is intended only for examples and software checks; it does not contain
#' participant data from a real study.
#'
#' @format A list with seven elements:
#' \describe{
#'   \item{A}{Binary group indicator.}
#'   \item{X}{Data frame with three baseline covariates.}
#'   \item{Z}{Data frame with four post-exposure covariates.}
#'   \item{M}{Continuous motion measure.}
#'   \item{Delta_M}{Binary motion-inclusion indicator.}
#'   \item{Delta_Y}{Binary outcome-observation and quality indicator.}
#'   \item{Y}{Numeric matrix with seven outcomes; the final column is a
#'     structural seed/self-correlation position.}
#' }
#' @usage data(data)
#' @keywords datasets
"data"
