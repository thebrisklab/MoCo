## Full MoCo analysis of the bundled example data
##
## This script reproduces the README tutorial: it fits MoCo to `moco_example`
## with three motion-density models (log-normal GLM, GAMLSS and HAL), fits
## four naive comparison methods, and compares all results with the simulated
## truth.
##
## Run time: with the default Super Learner library and five-fold cross-fitting,
## the GLM and GAMLSS fits each take about 15 to 16 minutes and the HAL fit
## about 20 minutes on one core (about 50 minutes in total); the naive methods
## take less than a second. It is not run during package checks; the examples
## in ?moco use a small, fast case instead.
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

## Settings shared by all fits: default Super Learner library, five-fold
## cross-fitting, and simultaneous tests at FWER 0.05
fit_moco <- function(...) {
  moco(
    X = d$X, Z = d$Z, A = d$A, M = d$M, Y = d$Y,
    Delta_M = d$Delta_M, Delta_Y = d$Delta_Y,
    cross_fit = TRUE, cv_folds = 5, seed_rng = 123,
    test = TRUE, fwer = 0.05, hypo_seed = 123, ...
  )
}

## 1. Log-normal GLM motion densities
fit_glm <- fit_moco(density_method = "glm")

## 2. GAMLSS: generalized gamma, penalized splines in age (X2) and IQ (Z2)
##    for the location, constant scale, shape nu estimated
fit_gamlss <- fit_moco(
  density_method = "gamlss",
  gamlss_options = list(mu = "smooth", smooth = c("X2", "Z2"), sigma = "1")
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
truth <- c(0, 0, 0, 0, -0.04849, -0.06805)
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

## Naive comparison methods. Naive LM: Y ~ A + X + M among Delta_Y = 1;
## removal LM: the same model after also removing Delta_M = 0; t tests: group
## mean differences in the same two samples. Each uses influence-function
## z-statistics and the same multiplier simultaneous test as moco().
naive_test <- function(effect, IF) {
  h <- hypo_test(
    list(adj_association = effect,
         eif_mat = lapply(seq_len(ncol(IF)), function(j) cbind(0, IF[, j])),
         cov_mat = lapply(seq_len(ncol(IF)), function(j) stats::cov(cbind(0, IF[, j])) / nrow(IF))),
    fwer = 0.05, method = "multiplier", n_sim = 10000L, seed = 123
  )
  list(est = effect, z = h$z_score, crit = h$fwer_critical_value)
}
naive_lm <- function(rows) {
  Yk <- d$Y[rows, 1:6]
  D <- model.matrix(~ A + X1 + X2 + X3 + M, data.frame(A = d$A, d$X, M = d$M)[rows, ])
  B <- qr.coef(qr(D), Yk); R <- Yk - D %*% B
  w <- drop(D %*% chol2inv(chol(crossprod(D) / nrow(D)))[, match("A", colnames(D))])
  naive_test(B["A", ], R * w)
}
naive_t <- function(rows) {
  Yk <- d$Y[rows, 1:6]; a <- d$A[rows]; p1 <- mean(a)
  m1 <- colMeans(Yk[a == 1, ]); m0 <- colMeans(Yk[a == 0, ])
  naive_test(m1 - m0, (a / p1) * sweep(Yk, 2, m1) - ((1 - a) / (1 - p1)) * sweep(Yk, 2, m0))
}
all_Y  <- d$Delta_Y == 1
remove <- d$Delta_Y == 1 & d$Delta_M == 1
naive <- list(naive_LM = naive_lm(all_Y), removal_LM = naive_lm(remove),
              naive_t = naive_t(all_Y), removal_t = naive_t(remove))

est_n <- sapply(naive, `[[`, "est"); z_n <- sapply(naive, `[[`, "z")
rownames(est_n) <- rownames(z_n) <- outcomes
cbind(truth = truth, round(est_n, 4))          # naive estimates
round(z_n, 2)                                  # z-statistics
sapply(naive, `[[`, "crit")                    # critical values
sweep(abs(z_n), 2, sapply(naive, `[[`, "crit"), ">")   # significant at FWER 0.05
round(sqrt(colMeans((est_n - truth)^2)), 3)    # RMSE over the six outcomes
