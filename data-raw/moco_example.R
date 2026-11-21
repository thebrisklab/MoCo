# -------------------
# MoCo package example data: data/moco_example.rda
# -------------------
# Simulated ASD vs non-ASD seed-based functional-connectivity study of n = 400
# children aged 8 to 14 years. The seed is the default mode network (DMN), and
# the outcomes are its Fisher-z connectivity with six Yeo networks.
#
# The data set contains:
#   X        baseline covariates: X1 sex (1 = male), X2 age (years), X3 handedness (1 = right)
#   A        diagnosis (1 = ASD, 0 = non-ASD)
#   Z        post-diagnosis covariates: Z1 ADOS (0 for non-ASD), Z2 full-scale IQ (FIQ),
#            Z3 stimulant medication, Z4 non-stimulant medication
#   M        head motion: mean framewise displacement (mm)
#   Delta_M  motion usability (1 = passed motion quality control)
#   Delta_Y  imaging-outcome availability (1 = available)
#   Y        400 x 7 Fisher-z connectivity matrix; rows with Delta_Y = 0 are NA, and the
#            seed-to-self column (Default) is all NA
#
# Generative model (expit(x) = 1 / (1 + exp(-x))):
#
#  (1) X1 ~ Bernoulli(0.75); X2 ~ Gamma(shape 25, rate 2.5), redrawn until it lies in
#      [8, 14]; X3 ~ Bernoulli(0.92).
#  (2) A ~ Bernoulli(expit(-0.11 + 0.71 X1 - 0.08 X2 - 0.19 X3)).
#  (3) Z1 = A * Poisson(11.86); Z2 ~ N(104.2, 17.4^2) if A = 1 and N(114.6, 11.6^2) if
#      A = 0; Z3 = A * Bernoulli(0.20); Z4 = A * Bernoulli(0.17).
#  (4) M ~ GG(mu, sigma, nu), the generalized gamma distribution of gamlss.dist, with
#      nu = -0.570 (a right tail heavier than the log-normal, nu = 0). With
#      g2 = (10 / X2)^2 - 1, g3 = (X2 / 10)^3 - 1 (fractional polynomial in age, powers
#      -2 and 3) and q = (Z2 - 110) / 15:
#        log mu    = -1.911 + 0.088 A + 0.60 g2 + 0.10 g3 + 0.10 X1 - 0.05 X3
#                    + 0.010 Z1 - 0.12 q + 0.05 q^2 - 0.10 Z3 - 0.05 Z4
#        log sigma = -0.203 - 0.096 A + 0.30 g2 + 0.05 g3 + 0.03 X1 - 0.03 X3
#                    + 0.005 Z1 - 0.05 q + 0.03 q^2 - 0.05 Z3 - 0.02 Z4
#      Younger children, children with higher ADOS and children with low IQ move more;
#      medicated children move slightly less. The four intercepts (log mu, log sigma, and
#      their A effects) and nu were calibrated on a synthetic population of 50,000 so
#      that the population 5th percentile, median and 95th percentile of M are about
#      0.052 / 0.174 / 0.90 mm (non-ASD) and 0.063 / 0.236 / 1.42 mm (ASD); see
#      calibrate_motion() below.
#  (5) P(Delta_M = 1) = 0.05 + 0.90 expit(2.658 - 2.8 log(M / 0.2) - 0.518 A
#                                         + 0.05 (X2 - 10.5)):
#      a smooth quality-control probability rather than a fixed cut-off, calibrated so
#      that about 80% of non-ASD and 67% of ASD children pass in the population (see
#      calibrate_qc() below). Delta_M is independent of Y given (M, A, X, Z).
#  (6) logit P(Delta_Y = 1) = 1.984 - 0.115 X1 - 0.0574 X2 + 0.0380 X3 - 0.0578 Z1
#                             + 0.0154 Z2 + 0.330 Z3.
#  (7) Y = outcome_mean() + E, with E = N(0, I_6) L (correlated errors, covariance
#      t(L) L). Motion enters the means through Mc = min(M, 1) and the hinge
#      h = max(Mc - 0.2, 0) with beta = 2, so the outcomes are nonlinear in motion and
#      the motion effect levels off above 1 mm. Only outcomes 5 (Limbic) and 6
#      (Frontoparietal) depend on A and Z.
#
# True MoCo group differences (A = 1 minus A = 0; motion standardized to
# p(m | Delta_M = 1, A = 0, X), Z marginalized given (A, M, X)), outcomes 1-7:
#   0, 0, 0, 0, -0.04849, -0.06805, NA
# Because M depends on Z, the truth is computed numerically by true_difference() below
# (about 10 minutes on one core; not run by default).
#
# Realized data (seed 123): 138 ASD; Delta_M = 1 for 78.6% of non-ASD and 65.9% of ASD
# children; Delta_Y = 1 for 372 children; median M 0.187 mm (non-ASD) and 0.221 mm (ASD).
#
# RNG: R default (Mersenne-Twister / Inversion / Rejection); data seed 123.
# Requires gamlss.dist (rGG, qGG, dGG). Run from the package root.
# -------------------

library(gamlss.dist)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# -------------------
# settings
# -------------------
n <- 400L          # sample size
seed <- 123L       # data seed
beta <- 2          # hinge coefficient of max(min(M, 1) - 0.2, 0) in the outcome means
M_cap <- 1         # motion above M_cap (mm) has the same effect on Y as M_cap

# motion model: intercept of log mu, intercept of log sigma, A effect on log mu,
# A effect on log sigma, and the GG shape nu (values from calibrate_motion())
motion_par <- c(log_mu0 = -1.911, log_sigma0 = -0.203, mu_A = 0.088, sigma_A = -0.096)
motion_nu <- -0.570
# quality-control model: intercept and A effect (values from calibrate_qc())
qc_par <- c(c0 = 2.658, cA = -0.518)
qc_floor <- 0.05; qc_ceiling <- 0.95

# Delta_Y model: intercept and coefficients of X1, X2, X3, Z1, Z2, Z3, Z4
b_DY <- c(1.98369441015941, -0.114997589466419, -0.0574406338144264, 0.0379817217336309,
          -0.057827471706101, 0.0154125907478968, 0.330170374268356, 0)

# loading matrix of the outcome errors (error covariance = crossprod(L))
L <- structure(c(-0.0579836993127972, 0.0923962425986226, -0.065496732195895, -0.0703380646194092, 0.0404827272535643, -0.0151935862169921,
  -0.159911224569557, -0.0913463631595361, -0.00781914709641775, -0.0351947520080434, -0.0150563668821247, -0.0404586661696976,
  -0.0956247682739167, 0.0693306744076057, -0.0341547928703659, 0.0313958560002284, -0.0768115088541821, 0.0133416665599556,
  -0.142720865659839, -0.0176185747949291, 0.0113502642976967, 0.0172020337382051, 0.0392264138498584, 0.0632554709482431,
  0.0420338597587603, -0.0559212927765671, -0.123316258021322, 0.0513111661107461, 0.0110113362154677, 0.0046560066317363,
  -0.0645533735127235, 0.0431279340904126, 0.0234035565095793, 0.0992356013624827, 0.0351623410760202, -0.0427116985290263), dim = c(6L, 6L))
# the error columns are permuted so that outcomes 5 and 6 carry the A effects
error_order <- c(1, 2, 3, 5, 4, 6)

outcome_names <- c("Visual", "Somatomotor", "DorsalAttention", "VentralAttention",
                   "Limbic", "Frontoparietal", "Default")

# -------------------
# model components
# -------------------
# X, A and Z
covariates <- function(n){
  X1 <- rbinom(n, 1, 0.75)                                  # sex (1 = male)
  X2 <- rgamma(n, 25, rate = 2.5)                           # age (years), in [8, 14]
  while(any(bad <- X2 < 8 | X2 > 14)) X2[bad] <- rgamma(sum(bad), 25, rate = 2.5)
  X3 <- rbinom(n, 1, 0.92)                                  # handedness (1 = right)
  A  <- rbinom(n, 1, plogis(-0.11 + 0.71*X1 - 0.08*X2 - 0.19*X3))   # ASD
  Z1 <- A*rpois(n, 11.86)                                   # ADOS (0 for non-ASD)
  Z2 <- ifelse(A == 1, rnorm(n, 104.2, 17.4), rnorm(n, 114.6, 11.6))  # FIQ
  Z3 <- A*rbinom(n, 1, 0.20)                                # stimulant medication
  Z4 <- A*rbinom(n, 1, 0.17)                                # non-stimulant medication
  data.frame(A, X1, X2, X3, Z1, Z2, Z3, Z4)
}

# log mu and log sigma of the generalized gamma motion model
motion_lp <- function(w, par = motion_par){
  g2 <- (10/w$X2)^2 - 1; g3 <- (w$X2/10)^3 - 1; q <- (w$Z2 - 110)/15
  with(w, list(
    log_mu    = par[["log_mu0"]] + par[["mu_A"]]*A + 0.60*g2 + 0.10*g3 + 0.10*X1 - 0.05*X3 +
                0.010*Z1 - 0.12*q + 0.05*q^2 - 0.10*Z3 - 0.05*Z4,
    log_sigma = par[["log_sigma0"]] + par[["sigma_A"]]*A + 0.30*g2 + 0.05*g3 + 0.03*X1 - 0.03*X3 +
                0.005*Z1 - 0.05*q + 0.03*q^2 - 0.05*Z3 - 0.02*Z4))
}

# probability of passing motion quality control
qc_prob <- function(M, A, X2, par = qc_par){
  qc_floor + (qc_ceiling - qc_floor)*plogis(par[["c0"]] - 2.8*log(M/0.2) + par[["cA"]]*A + 0.05*(X2 - 10.5))
}

# conditional means of Y (outcomes 1-6)
outcome_mean <- function(w){
  Mc <- pmin(w$M, M_cap); h <- pmax(Mc - 0.2, 0)
  with(w, cbind(
    -0.22 - 0.98*Mc - 0.06*X1 + 0.012*X2 + 0.03*X3 - beta*h,
    -0.20 + 0.92*Mc + 0.06*X1 - 0.009*X2 - 0.03*X3 + beta*h,
    -0.37 + 0.86*Mc + 0.04*X1 + 0.002*X2 + 0.04*X3 + beta*h,
     0.17 - 1.02*Mc - 0.06*X1 + 0.002*X2 + 0.04*X3 - beta*h,
    -0.20 - 0.03*A + 1.50*Mc - 0.61*Mc^2 + 0.02*X1 - 0.002*X2 + 0.03*X3 - 0.0005*Z1 + 0.0003*Z2 - 0.03*Z3 - 0.02*Z4 + beta*h,
    -0.16 - 0.05*A + 1.67*Mc - 0.64*Mc^2 + 0.03*X1 - 0.001*X2 + 0.02*X3 - 0.0005*Z1 + 0.0003*Z2 - 0.02*Z3 - 0.03*Z4 + beta*h))
}

# -------------------
# generate the example data
# -------------------
generate_moco_example <- function(n, seed){
  set.seed(seed)
  w <- covariates(n)
  # motion M (mean FD, mm)
  lp <- motion_lp(w)
  w$M <- rGG(n, mu = exp(lp$log_mu), sigma = exp(lp$log_sigma), nu = motion_nu)
  # Delta_Y (1 = usable imaging data)
  Delta_Y <- rbinom(n, 1, plogis(drop(as.matrix(cbind(1, w[, c("X1","X2","X3","Z1","Z2","Z3","Z4")])) %*% b_DY)))
  # Delta_M (1 = passed motion quality control)
  Delta_M <- rbinom(n, 1, qc_prob(w$M, w$A, w$X2))
  # Y: conditional means + correlated errors; 7th column is the seed-to-self position
  E <- (matrix(rnorm(n*6), n) %*% L)[, error_order]
  Y <- cbind(outcome_mean(w) + E, NA_real_)
  Y[Delta_Y == 0, ] <- NA_real_
  colnames(Y) <- outcome_names
  list(A = w$A, X = w[, c("X1","X2","X3")], Z = w[, c("Z1","Z2","Z3","Z4")], M = w$M,
       Delta_M = Delta_M, Delta_Y = Delta_Y, Y = Y)
}

moco_example <- generate_moco_example(n, seed)

# save the package data (uncomment to write)
# save(moco_example, file = "data/moco_example.rda", compress = "xz", version = 2)

# -------------------
# calibration of the motion and quality-control models (not run; reproduces
# motion_par, motion_nu and qc_par up to optimizer tolerance)
# -------------------
calibrate_motion <- function(N = 50000, seed = 20261001){
  set.seed(seed); P <- covariates(N); u <- runif(N)
  target <- log(c(0.052, 0.174, 0.90, 0.063, 0.236, 1.42))   # q05, median, q95 by group
  weight <- c(1, 10, 1, 1, 10, 1)
  quant <- function(th){
    par <- c(log_mu0 = th[1], log_sigma0 = th[2], mu_A = th[3], sigma_A = th[4])
    lp <- motion_lp(P, par); M <- qGG(u, exp(lp$log_mu), exp(lp$log_sigma), th[5])
    c(quantile(M[P$A == 0], c(.05, .5, .95)), quantile(M[P$A == 1], c(.05, .5, .95)))
  }
  loss <- function(th) if(th[5] < -3 || th[5] > 2) 1e6 else sum(weight*(log(quant(th)) - target)^2)
  optim(c(-1.75, -0.5, 0.2, 0.1, -1), loss, control = list(maxit = 1500, reltol = 1e-9))$par
}
calibrate_qc <- function(N = 50000, seed = 20261001){
  set.seed(seed); P <- covariates(N); lp <- motion_lp(P)
  P$M <- rGG(N, exp(lp$log_mu), exp(lp$log_sigma), motion_nu)
  loss <- function(p){ pr <- qc_prob(P$M, P$A, P$X2, c(c0 = p[1], cA = p[2]))
    (mean(pr[P$A == 0]) - 0.797)^2 + (mean(pr[P$A == 1]) - 0.672)^2 }
  optim(c(2.4, -0.4), loss)$par
}

# -------------------
# true MoCo group differences (not run; about 10 minutes)
# -------------------
# alpha_A + E_X int E[f(Z) | A = a, m, X] p(m | Delta_M = 1, A = 0, X) dm, a = 1 minus a = 0,
# where f(Z) is the Z part of the outcome-5 and outcome-6 means; the M and X parts cancel
# because the motion distribution is standardized to the same reference in both groups.
true_difference <- function(n_X = 1500, K = 1500, seed = 11){
  set.seed(seed); Xs <- covariates(n_X)
  m_grid <- seq(0.001, 2, length.out = 800)
  draw_Z <- function(a){ w <- covariates(3*K); w[w$A == a, ][1:K, ] }   # Z | A
  Z1s <- draw_Z(1); Z0s <- draw_Z(0)
  fZ <- function(z) cbind(-0.0005*z$Z1 + 0.0003*z$Z2 - 0.03*z$Z3 - 0.02*z$Z4,
                          -0.0005*z$Z1 + 0.0003*z$Z2 - 0.02*z$Z3 - 0.03*z$Z4)
  F1 <- fZ(Z1s); F0 <- fZ(Z0s)
  res <- t(sapply(seq_len(n_X), function(i){
    xi <- Xs[rep(i, K), c("X1","X2","X3")]
    l1 <- motion_lp(cbind(A = 1, xi, Z1s[, c("Z1","Z2","Z3","Z4")]))
    l0 <- motion_lp(cbind(A = 0, xi, Z0s[, c("Z1","Z2","Z3","Z4")]))
    D1 <- sapply(m_grid, function(m) dGG(rep(m, K), exp(l1$log_mu), exp(l1$log_sigma), motion_nu))
    D0 <- sapply(m_grid, function(m) dGG(rep(m, K), exp(l0$log_mu), exp(l0$log_sigma), motion_nu))
    ref <- colMeans(D0)*qc_prob(m_grid, 0, xi$X2[1]); ref <- ref/sum(ref)   # p(m | Delta_M = 1, A = 0, X)
    e1 <- crossprod(F1, D1)/matrix(pmax(colSums(D1), 1e-300), 2, length(m_grid), byrow = TRUE)
    e0 <- crossprod(F0, D0)/matrix(pmax(colSums(D0), 1e-300), 2, length(m_grid), byrow = TRUE)
    c(e1 %*% ref, e0 %*% ref)
  }))
  d <- colMeans(res)
  c(0, 0, 0, 0, c(-0.03, -0.05) + d[1:2] - d[3:4], NA)
}
