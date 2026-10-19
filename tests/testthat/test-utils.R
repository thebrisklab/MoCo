test_that("denominator densities use an absolute floor", {
  dens <- c(0, 1e-8, 1, 2, 3, NA)
  out <- MoCo:::.moco_floor_density(dens, 1e-8)
  expect_equal(out, c(1e-8, 1e-8, 1, 2, 3, NA))
})

test_that("density ratios can be capped explicitly", {
  ratio <- c(0, 2, 20, NA)
  expect_equal(
    MoCo:::.moco_cap_density_ratio(ratio, 10),
    c(0, 2, 10, NA)
  )
  expect_identical(MoCo:::.moco_cap_density_ratio(ratio, Inf), ratio)
})

test_that("EIF terms are zeroed by indicator and non-finite values are caught", {
  n <- 6
  A <- c(1, 1, 1, 0, 0, 0)
  Delta_Y <- c(1, 1, 0, 1, 0, 1)
  Delta_M <- c(1, 0, 1, 1, 1, 0)
  Y <- c(1, 2, NA, 3, NA, 4)
  args <- list(
    a = 1, A = A, Delta_Y = Delta_Y, Delta_M = Delta_M, Y = Y,
    gA = rep(0.5, n), gDM = rep(0.8, n), gDY_AXZ = rep(0.9, n),
    eta_AXZ = rep(1, n), eta_AXM = rep(1, n), xi_AX = rep(1, n),
    mu = rep(1, n), pMXD = rep(1, n), pMXZ = rep(1, n)
  )
  expect_true(all(is.finite(do.call(MoCo:::make_full_data_eif, args))))
  args$mu[1] <- NA
  expect_error(do.call(MoCo:::make_full_data_eif, args), "Non-finite")
})

test_that("learners from missing packages are dropped with a message", {
  lib <- c("SL.mean", "SL.glm")
  expect_identical(MoCo:::.moco_available_learners(lib), lib)
  skip_if(requireNamespace("xgboost", quietly = TRUE))
  expect_message(
    out <- MoCo:::.moco_available_learners(c("SL.mean", "SL.xgboost")),
    "Dropping"
  )
  expect_identical(out, "SL.mean")
})

test_that("the log-likelihood metalearner handles duplicate learners", {
  method <- MoCo:::tmp_method.CC_nloglik()
  Z <- cbind(rep(0.5, 10), rep(0.5, 10))
  fit <- suppressWarnings(method$computeCoef(
    Z = Z, Y = rep(c(0, 1), 5), libraryNames = c("a", "b"),
    obsWeights = rep(1, 10), control = list(trimLogit = 0.001),
    verbose = FALSE
  ))
  expect_length(fit$coef, 2)
  expect_equal(sum(fit$coef), 1)
})

test_that("support rules zero the target density where intended", {
  M <- c(0.10, 0.30, 0.08, NA)
  Delta_M <- c(1, 0, 0, 0)
  dens <- c(1, 2, 3, NA)
  expect_equal(
    MoCo:::.moco_apply_support_rule(dens, M, Delta_M, "delta_m", hal = FALSE),
    c(1, 0, 0, 0)
  )
  expect_identical(
    MoCo:::.moco_apply_support_rule(dens, M, Delta_M, "delta_m", hal = TRUE),
    dens
  )
  expect_equal(
    MoCo:::.moco_apply_support_rule(dens, M, Delta_M, "observed_range"),
    c(1, 0, 3, NA)
  )
  expect_equal(
    MoCo:::.moco_apply_support_rule(dens, M, Delta_M, "observed_range", support_upper = 0.05),
    c(0, 0, 0, NA)
  )
  expect_error(
    MoCo:::.moco_apply_support_rule(dens, M, rep(0, 4), "observed_range"),
    "Delta_M = 1"
  )
})
