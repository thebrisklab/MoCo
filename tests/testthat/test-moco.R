run_glm <- function(...) {
  e <- new.env()
  utils::data("data", package = "MoCo", envir = e)
  d <- e$data
  args <- list(
    X = d$X, Z = d$Z, A = d$A, M = d$M, Y = d$Y,
    Delta_M = d$Delta_M, Delta_Y = d$Delta_Y,
    SL_library = c("SL.mean", "SL.glm"),
    glm_formula = list(pMX = ".", pMXZ = "."),
    HAL_pMX = FALSE, HAL_pMXZ = FALSE,
    cross_fit = FALSE, test = TRUE, fwer = 0.05
  )
  extra <- list(...)
  args[names(extra)] <- extra
  suppressWarnings(do.call(moco, args))
}

example_data <- function() {
  e <- new.env()
  utils::data("data", package = "MoCo", envir = e)
  e$data
}

test_that("GLM run returns finite estimates with NA at the seed column", {
  fit <- run_glm()
  est <- as.matrix(fit$est)
  expect_equal(dim(est), c(2L, 7L))
  expect_true(all(is.na(est[, 7])))
  expect_true(all(is.finite(est[, -7])))
  expect_length(fit$adj_association, 7L)
  expect_true(is.na(fit$adj_association[7]))
  expect_true(all(is.finite(fit$adj_association[-7])))
  expect_true(all(is.finite(fit$z_score[-7])))
  expect_length(fit$fwer_critical_value, 1L)
})

test_that("moco() leaves the caller's random-number stream untouched", {
  set.seed(123)
  a <- stats::runif(1)
  set.seed(123)
  invisible(run_glm())
  b <- stats::runif(1)
  expect_identical(a, b)
})

test_that("thresh = 0.2 gives the same result as the supplied Delta_M", {
  d <- example_data()
  expect_identical(as.numeric(d$M < 0.2), as.numeric(d$Delta_M))
  fit_delta <- run_glm()
  fit_thresh <- run_glm(thresh = 0.2)
  expect_identical(fit_thresh$est, fit_delta$est)
  expect_identical(fit_thresh$adj_association, fit_delta$adj_association)
})

test_that("GAMLSS LOGNO agrees with the log-normal GLM densities", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  fit_glm <- run_glm()
  fit_gam <- run_glm(gamlss_options = list(
    family = "LOGNO", formulas = list(mu = "linear", sigma = "1")
  ))
  expect_identical(is.na(fit_gam$adj_association), is.na(fit_glm$adj_association))
  expect_lt(max(abs(fit_gam$adj_association - fit_glm$adj_association), na.rm = TRUE), 1e-3)
  expect_s3_class(fit_gam$gamlss_info, "data.frame")
})

test_that("default GG family fails on the example data with a clear message", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  expect_error(
    run_glm(gamlss_options = list(family = "GG")),
    "Try another family"
  )
})
