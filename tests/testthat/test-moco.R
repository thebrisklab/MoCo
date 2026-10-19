test_that("moco validates its inputs", {
  data("data", package = "MoCo", envir = environment())
  base_args <- list(
    X = data$X, Z = data$Z, A = data$A, M = data$M, Y = data$Y,
    Delta_M = data$Delta_M, Delta_Y = data$Delta_Y,
    pMX_method = "GLM", glm_formula = list(pMX = ".", pMXZ = ".")
  )
  expect_error(do.call(moco, modifyList(base_args, list(Delta_M = NULL))), "Delta_M or thresh")
  expect_error(do.call(moco, modifyList(base_args, list(cv_folds = 1))), "cv_folds")
  expect_error(do.call(moco, modifyList(base_args, list(density_floor = -1))), "density_floor")
  expect_error(do.call(moco, modifyList(base_args, list(density_floor = 0))), "density_floor")
  expect_error(do.call(moco, modifyList(base_args, list(density_ratio_cap = 0))), "density_ratio_cap")
  expect_error(do.call(moco, modifyList(base_args, list(A = data$A[-1]))), "same number")
  bad_M <- data$M
  bad_M[which(data$Delta_Y == 1)[1]] <- NA
  expect_error(do.call(moco, modifyList(base_args, list(M = bad_M))), "observed")
  # modifyList() drops empty lists, so set the element directly
  bad_formula <- base_args
  bad_formula$glm_formula <- list()
  expect_error(do.call(moco, bad_formula), "glm_formula")
})

test_that("the GLM path returns finite estimates with the right shape", {
  skip_on_cran()
  data("data", package = "MoCo", envir = environment())
  set.seed(99)
  before <- .Random.seed
  fit <- quiet_moco(
    X = data$X, Z = data$Z, A = data$A, M = data$M, Y = data$Y[, c(1, 7)],
    Delta_M = data$Delta_M, Delta_Y = data$Delta_Y,
    SL_library = c("SL.mean", "SL.glm"),
    pMX_method = "GLM", glm_formula = list(pMX = ".", pMXZ = "."),
    cv_folds = 2, test_n_sim = 1000
  )
  expect_identical(.Random.seed, before)
  expect_equal(dim(fit$est), c(2L, 2L))
  expect_true(all(is.finite(as.matrix(fit$est)[, 1])))
  expect_true(all(is.na(as.matrix(fit$est)[, 2])))
  expect_true(is.finite(fit$z_score[1]))
  expect_equal(fit$density_controls$density_floor, sqrt(.Machine$double.eps))
  expect_identical(fit$density_controls$density_ratio_cap, Inf)
})

test_that("support rules agree when Delta_M is a cutoff on M", {
  skip_on_cran()
  data("data", package = "MoCo", envir = environment())
  args <- list(
    X = data$X, Z = data$Z, A = data$A, M = data$M, Y = data$Y[, 1, drop = FALSE],
    Delta_M = data$Delta_M, Delta_Y = data$Delta_Y,
    SL_library = c("SL.mean", "SL.glm"),
    pMX_method = "GLM", glm_formula = list(pMX = ".", pMXZ = "."),
    cross_fit = FALSE, test = FALSE
  )
  fit_dm <- do.call(quiet_moco, args)
  fit_or <- do.call(quiet_moco, c(args, list(support_rule = "observed_range")))
  expect_identical(fit_dm$support_rule, "delta_m")
  expect_identical(fit_or$support_rule, "observed_range")
  expect_equal(fit_dm$est, fit_or$est)
  expect_error(do.call(moco, c(args, list(support_rule = "other"))))
})

test_that("repeated seeds average estimates and EIFs", {
  skip_on_cran()
  data("data", package = "MoCo", envir = environment())
  args <- list(
    X = data$X, Z = data$Z, A = data$A, M = data$M, Y = data$Y[, 1:2],
    Delta_M = data$Delta_M, Delta_Y = data$Delta_Y,
    SL_library = c("SL.mean", "SL.glm"),
    pMX_method = "GLM", glm_formula = list(pMX = ".", pMXZ = "."),
    cross_fit = TRUE, cv_folds = 2, test_n_sim = 2000
  )
  one <- do.call(quiet_moco, c(args, list(seed_rgn = 1)))
  two <- do.call(quiet_moco, c(args, list(seed_rgn = 2)))
  both <- do.call(quiet_moco, c(args, list(seed_rgn = 1:2)))
  expect_equal(as.matrix(both$est), (as.matrix(one$est) + as.matrix(two$est)) / 2)
  expect_true(all(is.finite(both$z_score)))
})

test_that("moco runs without motion or baseline covariates", {
  skip_on_cran()
  data("data", package = "MoCo", envir = environment())
  fit <- quiet_moco(
    X = NULL, Z = data$Z, A = data$A, M = NULL, Y = data$Y[, 1, drop = FALSE],
    Delta_M = data$Delta_M, Delta_Y = data$Delta_Y,
    SL_library = c("SL.mean", "SL.glm"),
    cross_fit = FALSE, test = FALSE
  )
  expect_identical(unname(fit$density_method), c("none", "none"))
  expect_true(all(is.finite(as.matrix(fit$est))))
  expect_error(
    moco(
      X = data$X, Z = data$Z, A = data$A, M = NULL, Y = data$Y[, 1, drop = FALSE],
      thresh = 0.2, Delta_Y = data$Delta_Y, test = FALSE
    ),
    "thresh requires M"
  )
})

test_that("GAMLSS densities with spline terms fit and predict", {
  skip_on_cran()
  data("data", package = "MoCo", envir = environment())
  # A soft, scrubbing-like inclusion indicator, so that the Delta_M = 1 motion
  # values are not hard-truncated.
  set.seed(7)
  Delta_M <- stats::rbinom(length(data$M), 1, stats::plogis(3 - 12 * data$M))
  fit <- quiet_moco(
    X = data$X, Z = data$Z, A = data$A, M = data$M, Y = data$Y[, 5:6],
    Delta_M = Delta_M, Delta_Y = data$Delta_Y,
    SL_library = c("SL.mean", "SL.glm"),
    pMX_method = "GAMLSS", gamlss_continuous_X = "X2",
    gamlss_continuous_Z = c("Z1", "Z2"),
    support_rule = "observed_range",
    cross_fit = FALSE, test = FALSE
  )
  expect_identical(unname(fit$density_method), c("GAMLSS", "GAMLSS"))
  expect_true(all(is.finite(as.matrix(fit$est))))
})
