# Efficient influence functions (EIF) of the MoCo estimand, used by one_step() and one_step_cross().

#' Evaluate EIF for stochastic interventional (in)direct effects
#'
#' @param a a Diagnosis group (0 or 1) for the group-specific mean.
#' @param A A binary vector of length n \code{number of participants} representing the diagnosis status.
#' @param Delta_Y A binary vector indicating whether the outcome \code{Y} is available and passes quality control.
#' @param Delta_M A binary vector indicating whether motion is available and meets the motion inclusion criteria.
#' @param Y A vector of continuous outcome of interest.
#' @param gA Estimate of P(A = a | X_i), i = 1,...,n.
#' @param gDM Estimate of P(Delta_M = 1 | A_i = 0, X_i), i = 1,...,n.
#' @param gDY_AXZ Estimate of P(Delta_Y = 1 | A_i, X_i, Z_i), i = 1,...,n.
#' @param eta_AXZ Estimate of the motion-standardized mean given (A_i, X_i, Z_i).
#' @param eta_AXM Estimate of the Z-standardized mean given (A_i, M_i, X_i).
#' @param xi_AX Estimate of E(eta_AXZ | A_i, X_i).
#' @param mu Estimate of E(Y | A_i, M_i, X_i, Z_i, Delta_Y = 1).
#' @param pMXD Estimate of p(M_i | A = 0, X_i, Delta_M = 1).
#' @param pMXZ Estimate of p(M_i | A_i, X_i, Z_i, Delta_Y = 1).
#'
#' @return An n-length vector of the estimated EIF evaluated on the observations.

make_full_data_eif <- function(a, A, Delta_Y, Delta_M, Y, gA, gDM, gDY_AXZ, eta_AXZ, eta_AXM, xi_AX, mu, pMXD, pMXZ){
  ipw_a <- as.numeric(A == a)/gA
  ipw_a_prime <- as.numeric(Delta_Y == 1&A == a)/(gDY_AXZ*gA)
  ipw_a_MD <- as.numeric(A == 0 & Delta_M == 1)/(((1 - gA)*I(a==1) + gA*I(a==0))*gDM)
  
  c_star <- pMXD/pMXZ
  
  p1 <- ipw_a_prime * c_star * (Y - mu)
  p2 <- ipw_a * (eta_AXZ - xi_AX)
  p3 <- ipw_a_MD * (eta_AXM - xi_AX)
  p4 <- xi_AX - mean(xi_AX)
  
  p1[is.na(p1)] <- 0
  p3[is.na(p3)] <- 0
  
  return(p1 + p2 + p3 + p4)
}

#' Compute the reduced EIF used by the missing-motion estimation branch.
#'
#' @param a Diagnosis group (0 or 1) for the group-specific mean.
#' @param A Observed diagnosis vector, aligned to the prediction rows.
#' @param Delta_Y Outcome-usability indicator for each row.
#' @param Y Observed outcome vector; unavailable outcomes may be missing.
#' @param gA Predicted P(A = a | X), evaluated for the supplied group a.
#' @param gDY_AXZ Predicted P(Delta_Y = 1 | A = a, X, Z).
#' @param xi_AX Baseline-covariate-standardized prediction at A = a.
#' @param mu Outcome prediction E(Y | A = a, M, X, Z, Delta_Y = 1).
#' @return An EIF vector with one value per participant.
#' @noRd

make_full_data_eif_easy <- function(a, A, Delta_Y, Y, gA, gDY_AXZ, xi_AX, mu) {
  ipw_a <- as.numeric(A == a) / gA
  ipw_a_prime <- as.numeric(Delta_Y == 1 & A == a) / (gDY_AXZ * gA)
  
  p1 <- ipw_a_prime * (Y - mu)
  p2 <- ipw_a * (mu - xi_AX)
  p3 <- xi_AX - mean(xi_AX)
  
  p1[is.na(p1)] <- 0
  
  return(p1 + p2 + p3)
}
