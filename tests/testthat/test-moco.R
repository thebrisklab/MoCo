test_that("moco() returns the documented components", {
  fit <- fast_moco()

  expect_named(fit, c("est", "adj_association", "z_score",
                      "fwer_critical_value", "significant_regions"))
  expect_s3_class(fit$est, "data.frame")
  expect_equal(dim(fit$est), c(2L, 2L))
  expect_equal(rownames(fit$est), c("est_A0", "est_A1"))
  expect_length(fit$adj_association, 2)
  expect_length(fit$z_score, 2)
  expect_length(fit$fwer_critical_value, 1)
  expect_type(fit$significant_regions, "logical")
})

test_that("the seed-to-self column is returned as NA", {
  fit <- fast_moco()

  expect_true(is.finite(fit$adj_association[1]))
  expect_true(is.na(fit$adj_association[2]))
  expect_true(is.na(fit$z_score[2]))
  expect_true(is.na(fit$significant_regions[2]))
  expect_true(all(is.na(fit$est[, 2])))
})

test_that("moco() outputs are internally consistent", {
  fit <- fast_moco()

  # adjusted association is the A = 1 mean minus the A = 0 mean
  # (fit$est has no column names, so it is indexed as a matrix)
  est <- as.matrix(fit$est)
  expect_equal(fit$adj_association[1], unname(est["est_A1", 1] - est["est_A0", 1]))
  # significance is |z| above the simultaneous critical value
  expect_equal(fit$significant_regions[1],
               abs(fit$z_score[1]) > fit$fwer_critical_value)
  # one outcome: the simultaneous critical value is close to the normal quantile
  expect_equal(fit$fwer_critical_value, stats::qnorm(0.975), tolerance = 0.05)
})

test_that("moco() is in the neighbourhood of the simulated truth", {
  # True Frontoparietal difference is -0.06805 (data-raw/moco_example.R). This
  # small-library fit is not expected to match it closely; the check only
  # catches gross errors such as a flipped sign or a wrong scale.
  fit <- fast_moco()
  expect_lt(abs(fit$adj_association[1] - (-0.06805)), 0.1)
})

test_that("moco() is reproducible and leaves the caller's RNG alone", {
  set.seed(42)
  before <- .Random.seed
  fit1 <- fast_moco()
  expect_identical(.Random.seed, before)

  fit2 <- fast_moco()
  expect_identical(fit1, fit2)
})

test_that("test = FALSE returns estimates only", {
  fit <- fast_moco(test = FALSE)
  expect_named(fit, c("est", "adj_association"))
})

test_that("moco() runs without cross-fitting and with a vector outcome", {
  d <- moco_example
  fit <- fast_moco(cross_fit = FALSE, Y = d$Y[, 6])
  expect_length(fit$adj_association, 1)
  expect_true(is.finite(fit$adj_association))
  expect_true(is.finite(fit$z_score))
})

test_that("several fwer levels give one critical value and one row each", {
  fit <- fast_moco(fwer = c(0.01, 0.05, 0.10))
  expect_length(fit$fwer_critical_value, 3)
  expect_true(all(diff(fit$fwer_critical_value) < 0))
  expect_equal(dim(fit$significant_regions), c(3L, 2L))
})

test_that("thresh defines Delta_M when Delta_M is NULL", {
  d <- moco_example
  # M and Delta_Y must agree with the threshold rule on missingness only, so
  # the example data can be reused with a threshold-based Delta_M
  fit <- fast_moco(Delta_M = NULL, thresh = 0.3)
  expect_true(is.finite(fit$adj_association[1]))

  same <- fast_moco(Delta_M = as.numeric(d$M <= 0.3), thresh = 0.3)
  expect_identical(fit, same)
})

test_that("moco() rejects invalid inputs with informative errors", {
  d <- moco_example
  A_na <- d$A; A_na[1] <- NA
  DM_bad <- d$Delta_M; DM_bad[1] <- 2
  DY_na <- d$Delta_Y; DY_na[1] <- NA

  expect_error(fast_moco(A = A_na), "missing values in A, X, or Z")
  expect_error(fast_moco(Delta_M = NULL), "Delta_M or thresh is required")
  expect_error(fast_moco(Delta_M = DM_bad), "0/1 vector")
  expect_error(fast_moco(Delta_M = d$Delta_M[-1]), "length n")
  expect_error(fast_moco(Delta_Y = DY_na), "0/1 vector")
  expect_error(fast_moco(thresh = -1), "one positive number")
  expect_error(fast_moco(thresh = 0.3), "not equal to 1\\(M <= thresh\\)")
  # GitHub issue #1: M = NULL, Delta_M = NULL, thresh given
  expect_error(fast_moco(M = NULL, Delta_M = NULL, thresh = 0.3),
               "thresh requires M")
})

test_that("with M = NULL, Delta_Y alone defines the missing data", {
  fit <- fast_moco(M = NULL, Delta_M = NULL)
  expect_true(all(is.finite(fit$est[, 1])))
  expect_true(is.finite(fit$z_score[1]))

  # no motion densities are fitted, so the density method does not matter
  expect_identical(fit, fast_moco(M = NULL, Delta_M = NULL, density_method = "gamlss"))

  # without cross-fitting
  nocf <- fast_moco(M = NULL, Delta_M = NULL, cross_fit = FALSE)
  expect_true(is.finite(nocf$adj_association[1]))

  # Delta_M is not accepted without M
  expect_error(fast_moco(M = NULL), "Delta_M is not used when M = NULL")
})
