# -------------------
# data simulation: MoCo package example data (moco_example.rda)
# -------------------
# Generates the bundled MoCo example data `moco_example` (data/moco_example.rda): a
# simulated ASD vs non-ASD seed-based functional-connectivity study of n = 400 children
# motivated by ABIDE. The design follows the package example data of Ran et al. (2026),
# Supplement Section 9.2 (Example_Data/Code/sim_package_data.R). Only the motion M and the
# motion-QC indicator Delta_M are new, so the GLM, GAMLSS and HAL motion-density options of
# moco() can be compared.
#
# RNG: R default (Mersenne-Twister / Inversion / Rejection).
# seed = 1004 (data); calibration_seed = 20261001 (calibration population, n = 50,000).
# Requires gamlss.dist (qGG, rGG). Run from the package root.
#
# Same as the original package example (sim_package_data.R):
#  (1) X, A and Z: distributions and coefficients.
#  (2) Delta_Y: the logistic model in X1, X2, X3, Z1-Z4 (b_DY; the Z4 coefficient is set
#      to exactly 0, it was 8e-17).
#  (3) Y: the six conditional means, the error loading matrix L, the column order (sim.R
#      columns 1, 2, 3, 5, 4, 6, so the A effects are in outcomes 5 and 6), Y[Delta_Y == 0, ]
#      = NA, and a 7th all-NA seed-to-self column (DMN with itself).
#  (4) n = 400 and the true group differences.
#
# Changes from the original package example:
#  (1) Motion M (mean FD, mm) is generalized gamma, M ~ GG(mu, sigma, nu = 1)
#      (gamlss.dist parameterization), instead of log-normal; mu and sigma depend on all
#      X and Z. With a = (X2 - 10.5)/2 and q = (Z2 - 110)/15:
#        log mu    = b0 + 0.15 A + 0.08 X1 - 0.05 X3 - 0.25 a + 0.12 a^2 + 0.12 Z1/12
#                       - 0.15 tanh(q) - 0.10 Z3 - 0.05 Z4
#        log sigma = g0 + 0.30 A + 0.05 X1 - 0.05 X3 - 0.15 a + 0.15 Z1/12
#                       - 0.05 tanh(q) - 0.10 Z3 - 0.05 Z4
#      More, and more variable, motion in ASD, at higher ADOS and at younger ages, and less
#      at higher IQ (Nebel et al. 2022, NeuroImage); medicated children move slightly less
#      (illustrative effects).
#      g0 and b0 are calibrated so that the population q95 / median of M is 0.55 / 0.18 mm.
#  (2) Delta_M (1 = usable motion) is probabilistic instead of a cut (M < 0.2):
#        P(Delta_M = 1) = 0.05 + 0.90 expit(c0 - 2.8 log(M / 0.2) - 0.20 A + 0.10 a),
#      with c0 calibrated to 80% usable in the population (about 84% non-ASD and 73% ASD;
#      84% TD and 72% ASD under lenient QC in Nebel et al. 2022). Delta_M is independent of
#      Y given (M, A, X, Z), as MoCo assumes.
#  (3) The random-number stream differs (seed 1004; M is drawn before Delta_Y), so the
#      realized X, A, Z and Y are not those of the original data.rda.
#
# Realized data (seed = 1004): 153 ASD; Delta_M = 1 for 79% (86% non-ASD, 68% ASD);
# Delta_Y = 1 for 93%; median M 0.16 mm (non-ASD) and 0.22 mm (ASD).
#
# True group differences (A = 1 minus A = 0) for outcomes 1-7:
#   0, 0, 0, 0, -0.04845, -0.06815, NA
#   (the outcome means have no A x M or Z x M terms, so the truth does not depend on the
#   motion model or on Delta_M.)
# -------------------

library(gamlss.dist)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# -------------------
# parameters from the original package example (sim_package_data.R)
# -------------------
# loading matrix of the outcome errors (mat_cov = crossprod(L)), in sim.R column order
L <- structure(c(-0.0579836993127972, 0.0923962425986226, -0.065496732195895, -0.0703380646194092, 0.0404827272535643, -0.0151935862169921,
  -0.159911224569557, -0.0913463631595361, -0.00781914709641775, -0.0351947520080434, -0.0150563668821247, -0.0404586661696976,
  -0.0956247682739167, 0.0693306744076057, -0.0341547928703659, 0.0313958560002284, -0.0768115088541821, 0.0133416665599556,
  -0.142720865659839, -0.0176185747949291, 0.0113502642976967, 0.0172020337382051, 0.0392264138498584, 0.0632554709482431,
  0.0420338597587603, -0.0559212927765671, -0.123316258021322, 0.0513111661107461, 0.0110113362154677, 0.0046560066317363,
  -0.0645533735127235, 0.0431279340904126, 0.0234035565095793, 0.0992356013624827, 0.0351623410760202, -0.0427116985290263), dim = c(6L, 6L))
# coefficients of the Delta_Y model:
# logit P(Delta_Y = 1) = b0 + b1*X1 + b2*X2 + b3*X3 + b4*Z1 + b5*Z2 + b6*Z3 + b7*Z4
b_DY <- c(1.98369441015941, -0.114997589466419, -0.0574406338144264, 0.0379817217336309, -0.057827471706101, 0.0154125907478968,
          0.330170374268356, 0)

# -------------------
# design settings
# -------------------
design <- list(nu = 1,                 # GG shape of the motion distribution
               target_median = 0.18,   # population median of M (mm)
               target_q95 = 0.55,      # population 95th percentile of M (mm)
               target_QC = 0.80,       # population proportion with Delta_M = 1
               QC_floor = 0.05,        # lower bound of P(Delta_M = 1)
               QC_ceiling = 0.95,      # upper bound of P(Delta_M = 1)
               calibration_seed = 20261001L, calibration_n = 50000L,
               n = 400L, seed = 1004L)

# -------------------
# model components
# -------------------
# X (demographic confounders), A (diagnosis) and Z (behavioural phenotypes), as in the original example
covariates <- function(n){
  # sex: 1 = male, 0 = female
  X1 <- rbinom(n, 1, 0.75)
  # age: Gamma(shape 25, rate 2.5), redrawn until all values lie in [8, 14]
  X2 <- rgamma(n, 25, rate = 2.5)
  while(any(bad <- X2 < 8 | X2 > 14)) X2[bad] <- rgamma(sum(bad), 25, rate = 2.5)
  # handedness: 1 = right-handed, 0 = left-handed
  X3 <- rbinom(n, 1, 0.92)
  # disease status: 1 = ASD, 0 = non-ASD
  A <- rbinom(n, 1, plogis(-0.11 + 0.71*X1 - 0.08*X2 - 0.19*X3))
  # ADOS score (0 for non-ASD)
  Z1 <- A*rpois(n, 11.86)
  # FIQ score
  Z2 <- ifelse(A == 1, rnorm(n, 104.2, 17.4), rnorm(n, 114.6, 11.6))
  # stimulant and non-stimulant medication (0 for non-ASD)
  Z3 <- A*rbinom(n, 1, 0.20)
  Z4 <- A*rbinom(n, 1, 0.17)
  data.frame(A, X1, X2, X3, Z1, Z2, Z3, Z4)
}
# linear predictors of log mu and log sigma of M (without intercepts) and of the QC model
predictors <- function(w){
  a <- (w$X2 - 10.5)/2; q <- (w$Z2 - 110)/15
  list(mu = 0.15*w$A + 0.08*w$X1 - 0.05*w$X3 - 0.25*a + 0.12*a^2 + 0.12*w$Z1/12 - 0.15*tanh(q) - 0.10*w$Z3 - 0.05*w$Z4,
       sigma = 0.30*w$A + 0.05*w$X1 - 0.05*w$X3 - 0.15*a + 0.15*w$Z1/12 - 0.05*tanh(q) - 0.10*w$Z3 - 0.05*w$Z4,
       QC = if(!is.null(w$M)) -2.8*log(w$M/0.2) - 0.20*w$A + 0.10*a)
}
# bounded logistic probability of usable motion
qc_prob <- function(eta) design$QC_floor + (design$QC_ceiling - design$QC_floor)*plogis(eta)
# conditional means of Y, in the final column order (outcomes 1-6)
outcome_mean <- function(w) with(w, cbind(
  -0.22 - 0.98*M - 0.06*X1 + 0.012*X2 + 0.03*X3,
  -0.20 + 0.92*M + 0.06*X1 - 0.009*X2 - 0.03*X3,
  -0.37 + 0.86*M + 0.04*X1 + 0.002*X2 + 0.04*X3,
   0.17 - 1.02*M - 0.06*X1 + 0.002*X2 + 0.04*X3,
  -0.20 - 0.03*A + 1.50*M - 0.61*M^2 + 0.02*X1 - 0.002*X2 + 0.03*X3 - 0.0005*Z1 + 0.0003*Z2 - 0.03*Z3 - 0.02*Z4,
  -0.16 - 0.05*A + 1.67*M - 0.64*M^2 + 0.03*X1 - 0.001*X2 + 0.02*X3 - 0.0005*Z1 + 0.0003*Z2 - 0.02*Z3 - 0.03*Z4))

# -------------------
# calibrate the intercepts on an independent synthetic population (before any method is fitted)
# -------------------
set.seed(design$calibration_seed)
cal <- covariates(design$calibration_n); lp <- predictors(cal); u <- runif(nrow(cal))
# g0 (log sigma intercept): match the q95 / median ratio of M
ratio <- function(g0){ m <- qGG(u, mu = exp(lp$mu), sigma = exp(g0 + lp$sigma), nu = design$nu)
  unname(quantile(m, 0.95)/median(m)) - design$target_q95/design$target_median }
design$log_sigma0 <- uniroot(ratio, c(log(0.2), log(2)), tol = 1e-7)$root
# b0 (log mu intercept): match the median of M
m0 <- qGG(u, mu = exp(lp$mu), sigma = exp(design$log_sigma0 + lp$sigma), nu = design$nu)
design$log_mu0 <- log(design$target_median/median(m0)); cal$M <- m0*exp(design$log_mu0)
# c0 (QC intercept): match the proportion with usable motion
qclp <- predictors(cal)$QC
design$QC_intercept <- uniroot(function(c0) mean(qc_prob(c0 + qclp)) - design$target_QC, c(-15, 15))$root

# -------------------
# generate the example data
# -------------------
gen <- function(n, seed){
  set.seed(seed)
  # X, A, Z
  w <- covariates(n); lp <- predictors(w)
  # motion M (mean FD, mm): generalized gamma
  w$M <- rGG(n, mu = exp(design$log_mu0 + lp$mu), sigma = exp(design$log_sigma0 + lp$sigma), nu = design$nu)
  # Delta_Y (1 = usable imaging data)
  Delta_Y <- rbinom(n, 1, plogis(drop(as.matrix(cbind(1, w[, c("X1","X2","X3","Z1","Z2","Z3","Z4")])) %*% b_DY)))
  # Delta_M (1 = usable motion): probabilistic QC
  Delta_M <- rbinom(n, 1, qc_prob(design$QC_intercept + predictors(w)$QC))
  # Y: conditional means + correlated errors (errors reordered to the final column order)
  E <- (matrix(rnorm(n*6), n) %*% L)[, c(1, 2, 3, 5, 4, 6)]
  Y <- cbind(outcome_mean(w) + E, NA_real_)
  Y[Delta_Y == 0, ] <- NA_real_
  colnames(Y) <- c("Visual", "Somatomotor", "DorsalAttention", "VentralAttention", "Limbic", "Frontoparietal", "Default")
  list(A = w$A, X = w[, c("X1","X2","X3")], Z = w[, c("Z1","Z2","Z3","Z4")], M = w$M,
       Delta_M = Delta_M, Delta_Y = Delta_Y, Y = Y)
}

# package data object ----------------------
moco_example <- gen(design$n, design$seed)

# save simulated data (uncomment to write)
# save(moco_example, file = "data/moco_example.rda", compress = "xz", version = 2)
