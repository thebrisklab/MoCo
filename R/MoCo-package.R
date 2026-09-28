#' MoCo: motion-controlled brain-phenotype estimation
#'
#' MoCo provides nonparametric one-step estimators for group-specific means and
#' contrasts under motion-related selection, with optional cross-fitting and
#' simultaneous family-wise error-rate control.
#'
#' @section Main functions:
#' - [moco()] fits MoCo using HAL, log-normal GLM, or generalized-gamma
#'   GAMLSS motion densities.
#' - [hypo_test()] applies EIF-based simultaneous confidence-band tests.
#'
#' @keywords internal
#' @importFrom stats binomial cov dnorm gaussian sd
"_PACKAGE"

utils::globalVariables(c("RS", "CG", "mixed"))
