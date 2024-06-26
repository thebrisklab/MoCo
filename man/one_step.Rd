% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/one_step.R
\name{one_step}
\alias{one_step}
\title{Non-parametric efficient motion-controlled brain phenotype estimators}
\usage{
one_step(
  X,
  Z,
  A,
  M,
  Y,
  Delta_M,
  thresh = NULL,
  Delta_Y,
  SL_library = c("SL.earth", "SL.glmnet", "SL.gam", "SL.glm", "SL.glm.interaction",
    "SL.step", "SL.step.interaction", "SL.xgboost", "SL.ranger", "SL.mean"),
  SL_library_customize = list(gA = NULL, gDM = NULL, gDY_AX = NULL, gDY_AXZ = NULL,
    mu_AMXZ = NULL, eta_AXZ = NULL, eta_AXM = NULL, xi_AX = NULL),
  glm_formula = list(gA = NULL, gDM = NULL, gDY_AX = NULL, gDY_AXZ = NULL, mu_AMXZ =
    NULL, eta_AXZ = NULL, eta_AXM = NULL, xi_AX = NULL, pMX = NULL, pMXZ = NULL),
  HAL_pMX = TRUE,
  HAL_pMXZ = TRUE,
  HAL_options = list(max_degree = 3, lambda_seq = exp(seq(-1, -10, length = 100)),
    num_knots = c(1000, 500, 250)),
  seed = 1,
  ...
)
}
\arguments{
\item{X}{A dataframe or matrix containing demographic confounders that would ideally be balanced in a randomized controlled trial.}

\item{Z}{A dataframe or matrix of covariates representing brain phenotypes.}

\item{A}{A binary vector of length n (number of participants), serving as a group indicator, such as diagnosis group or control group.}

\item{M}{A numeric vector of length n representing continuous motion values for each participant.}

\item{Y}{A matrix of dimension n x p, where n is the number of participants, and p is the number of regions of interest.
If it represents seed-based association measures:
Each (i, j) element denotes participant i's association measure between the seed region and region j.
The column representing the association measure of the seed region with itself should be filled with NA values to indicate its position.
If it represents other types of association measures:
Each (i, j) element denotes participant i's association measure between two brain regions of interest, such as the upper diagonal part of the functional connectivity matrix.
No NA values are allowed in \code{Y} in this case.}

\item{Delta_M}{A binary vector of length n indicating whether motion is available and meets inclusion criteria.
If motion meets inclusion criteria for analysis, set \code{Delta_M} = 1; otherwise, set \code{Delta_M} = 0.}

\item{thresh}{A numeric value used to threshold M to produce \code{Delta_M}. One can specify either \code{Delta_M} or thresh.}

\item{Delta_Y}{A binary vector indicating the non-missingness and whether the brain image data \code{Y} passes quality control after preprocessing.
Set \code{Delta_Y = 1} if \code{Y} is usable; otherwise, set \code{Delta_Y = 0}.}

\item{SL_library}{SuperLearner library for estimating nuisance regressions.
Defaults to c("SL.earth","SL.glmnet","SL.gam","SL.glm", "SL.glm.interaction", "SL.step","SL.step.interaction","SL.xgboost","SL.ranger","SL.mean") if not specified.}

\item{SL_library_customize}{Customize SuperLearner library for estimating each nuisance regression.
- \code{gA}: SuperLearner library for estimating the propensity score.
- \code{gDM}: SuperLearner library for estimating the probability P(Delta_M = 1 | A, X).
- \code{gDY_AX}: SuperLearner library for estimating the probability P(Delta_Y = 1 | A, X).
- \code{gDY_AXZ}: SuperLearner library for estimating the probability P(Delta_Y = 1 | A, X, Z).
- \code{mu_AMXZ}: SuperLearner library for estimating the outcome regression E(Y | Delta_Y = 1, A, M, X, Z).
- \code{eta_AXZ}: SuperLearner library for estimating E(mu_AMXZ pMXD / pMXZD | A, X, Z, Delta_M = 1).
- \code{eta_AXM}: SuperLearner library for estimating E(mu_AMXZ pMX/pMXZ gDY_AX/gDY_AXZ | A, M, X, Delta_Y = 1).
- \code{xi_AX}: SuperLearner library for estimating E(eta_AXZ | A, X).}

\item{glm_formula}{All glm formulas default to NULL, indicating SuperLearner will be used for nuisance regressions.
- \code{gA}: GLM formula for estimating the propensity score.
- \code{gDM}: GLM formula for estimating the probability P(Delta_M = 1 | A, X).
- \code{gDY_AX}: GLM formula for estimating the probability P(Delta_Y = 1 | A, X).
- \code{gDY_AXZ}: GLM formula for estimating the probability P(Delta_Y = 1 | A, X, Z).
- \code{mu_AMXZ}: GLM formula for estimating the outcome regression E(Y | Delta_Y = 1, A, M, X, Z).
- \code{eta_AXZ}: GLM formula for estimating E(mu_AMXZ pMXD / pMXZD | A, X, Z, Delta_M = 1).
- \code{eta_AXM}: GLM formula for estimating E(mu_AMXZ pMX/pMXZ gDY_AX/gDY_AXZ | A, M, X, Delta_Y = 1).
- \code{xi_AX}: GLM formula for estimating E(eta_AXZ | A, X).
- \code{pMX}: GLM formula for estimating p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M = 1), assuming M follows a log normal distribution.
- \code{pMXZ}: GLM formula for estimating p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M = 1), assuming M follows a log normal distribution.}

\item{HAL_pMX}{Specifies whether to estimate p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M=1) using the highly adaptive lasso conditional density estimation method.
Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMX} option in \code{glm_formula}, such as \code{pMX = "."}.}

\item{HAL_pMXZ}{Specifies whether to estimate p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M=1) using the highly adaptive lasso conditional density estimation method.
Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMXZ} option in \code{glm_formula}, such as \code{pMXZ = "."}.}

\item{HAL_options}{Additional options for highly adaptive lasso (HAL) method.
- \code{max_degree}: The highest order of interaction terms for generating basis functions (passed to \code{haldensify}).
- \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression (passed to \code{haldensify}).
- \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions (passed to \code{haldensify}).}
}
\value{
A list with named entries
\describe{
\item{est}{A two times p matrix showing the one-step estimators of the control group and the disease group for each functional connectivity of interest, respectively.}
\item{adj_association}{A p-length vector showing the motion-controlled association for each functional connectivity of interest, respectively.}
\item{eif_mat}{A p-length list of the estimated EIF evaluated on the observations of the control and disease group, respectively.}
\item{cov_mat}{A p-length list of the estimated covariance matrix of the one-step estimators.}
}
}
\description{
Non-parametric efficient motion-controlled brain phenotype estimators
}
\details{
Compute non-parametric efficient motion-controlled brain phenotype estimators (MoCo) using one-step estimation
}
