## Full MoCo analysis of the bundled example data
##
## This script reproduces the README tutorial: it fits MoCo to `moco_example`
## with three motion-density models (log-normal GLM, GAMLSS and HAL) and
## compares the results with the simulated truth.
##
## Run time: with the default Super Learner library and five-fold cross-fitting,
## the GLM and GAMLSS fits each take about 14 minutes and the HAL fit about
## 21 minutes on one core (about 50 minutes in total). It is not run during
## package checks; the examples in ?moco use a small, fast case instead.
##
## Open this file with
##   file.show(system.file("examples", "full_analysis.R", package = "MoCo"))
## or run it with
##   source(system.file("examples", "full_analysis.R", package = "MoCo"), echo = TRUE)
##
## Results may vary slightly with R, dependency and optimizer versions.

library(MoCo)
data(moco_example)
d <- moco_example

## Settings shared by the three fits: default Super Learner library,
## five-fold cross-fitting and simultaneous tests at FWER 0.05
fit_moco <- function(...) {
  moco(
    X = d$X, Z = d$Z, A = d$A, M = d$M, Y = d$Y,
    Delta_M = d$Delta_M, Delta_Y = d$Delta_Y,
    cross_fit = TRUE, cv_folds = 5, seed_rng = 1,
    test = TRUE, fwer = 0.05, ...
  )
}

## 1. Log-normal GLM motion densities
fit_glm <- fit_moco(density_method = "glm")

## 2. GAMLSS: generalized gamma, penalized splines in age (X2) and IQ (Z2)
##    for the location, scale linear in all covariates, shape nu estimated
fit_gamlss <- fit_moco(
  density_method = "gamlss",
  gamlss_options = list(mu = "smooth", smooth = c("X2", "Z2"), sigma = "linear")
)

## 3. HAL with a lighter basis than the default HAL_options
fit_hal <- fit_moco(
  density_method = "hal",
  HAL_options = list(
    max_degree = 2,
    lambda_seq = exp(seq(-1, -10, length = 40))[1:38],
    num_knots = c(30, 10, 5)
  )
)

## Comparison with the simulated truth (data-raw/moco_example.R).
## Column 7 of Y is the seed-to-self position and is NA by construction.
fits  <- list(GLM = fit_glm, GAMLSS = fit_gamlss, HAL = fit_hal)
truth <- c(0, 0, 0, 0, -0.04845, -0.06815)
outcomes <- c("Visual", "Somatomotor", "DorsalAttention", "VentralAttention",
              "Limbic", "Frontoparietal")

est <- sapply(fits, function(f) f$adj_association[1:6])
z   <- sapply(fits, function(f) f$z_score[1:6])
sig <- sapply(fits, function(f) f$significant_regions[1:6])
rownames(est) <- rownames(z) <- rownames(sig) <- outcomes

cbind(truth = truth, round(est, 4))  # adjusted associations
round(z, 2)                          # z-statistics
sig                                  # significant at FWER 0.05
sapply(fits, function(f) f$fwer_critical_value)
round(sqrt(colMeans((est - truth)^2)), 3)  # RMSE over the six outcomes

## Density-ratio weights pMXD / pMXZ: large weights or a small
## effective-sample-size ratio indicate a poorly fitting motion density
lapply(fits, function(f) f$weight_diagnostics)

## GAMLSS fit of each density in each fold (algorithm, convergence, nu, fallback)
fit_gamlss$density_log
