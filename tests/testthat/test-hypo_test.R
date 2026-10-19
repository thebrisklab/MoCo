test_that("hypo_test returns correctly shaped output", {
  res <- hypo_test(toy_result(), fwer = c(0.05, 0.1), n_sim = 2000)
  expect_named(res, c("z_score", "conf_band", "significant_regions", "correlation_matrix"))
  expect_length(res$z_score, 3)
  expect_equal(dim(res$significant_regions), c(2L, 3L))
  expect_true(res$conf_band[1] > res$conf_band[2])
})

test_that("one outcome gives a two-sided normal critical value", {
  res <- hypo_test(toy_result(p = 1), fwer = 0.05, n_sim = 20000)
  expect_equal(unname(res$conf_band), stats::qnorm(0.975), tolerance = 0.05)
})

test_that("hypo_test is reproducible and restores the caller's RNG state", {
  tr <- toy_result()
  set.seed(42)
  before <- .Random.seed
  r1 <- hypo_test(tr, n_sim = 1000, seed = 7)
  expect_identical(.Random.seed, before)
  r2 <- hypo_test(tr, n_sim = 1000, seed = 7)
  expect_identical(r1$conf_band, r2$conf_band)
})

test_that("hypo_test validates its inputs", {
  expect_error(hypo_test(toy_result(), fwer = 1.2), "fwer")
  expect_error(hypo_test(toy_result(), n_sim = 0), "n_sim")
  expect_error(hypo_test(list(adj_association = 1)), "missing")
})
