toy_result <- function(n = 200, p = 3, seed = 1) {
  set.seed(seed)
  eif <- lapply(seq_len(p), function(j) {
    a <- rnorm(n)
    cbind(a, a + rnorm(n))
  })
  list(adj_association = seq(0.05, 0.3, length.out = p),
       eif_mat = eif,
       cov_mat = lapply(eif, function(e) stats::cov(e) / n))
}

test_that("z-scores are the estimate over the contrast standard error", {
  toy <- toy_result()
  h <- hypo_test(toy, n_sim = 2000)

  se <- vapply(toy$eif_mat, function(e) stats::sd(e[, 2] - e[, 1]) / sqrt(nrow(e)),
               numeric(1))
  expect_equal(h$z_score, toy$adj_association / se)
})

test_that("hypo_test() returns one critical value and one row per fwer level", {
  h <- hypo_test(toy_result(), fwer = c(0.01, 0.05, 0.10), n_sim = 2000)

  expect_length(h$fwer_critical_value, 3)
  expect_true(all(diff(h$fwer_critical_value) < 0))
  expect_equal(dim(h$significant_regions), c(3L, 3L))
  expect_equal(h$significant_regions[2, ], abs(h$z_score) > h$fwer_critical_value[2])
})

test_that("critical values behave as simultaneous quantiles", {
  # one outcome: the normal quantile
  one <- hypo_test(toy_result(p = 1), n_sim = 20000)
  expect_equal(one$fwer_critical_value, stats::qnorm(0.975), tolerance = 0.03)

  # three independent outcomes: between the unadjusted and Bonferroni quantiles
  three <- hypo_test(toy_result(n = 2000, p = 3), n_sim = 20000)
  expect_gt(three$fwer_critical_value, stats::qnorm(0.975))
  expect_lt(three$fwer_critical_value, stats::qnorm(1 - 0.05 / 6) + 0.05)
})

test_that("the multiplier and mvn methods agree up to Monte Carlo error", {
  toy <- toy_result()
  mult <- hypo_test(toy, method = "multiplier", n_sim = 20000)
  mvn  <- hypo_test(toy, method = "mvn", n_sim = 20000)

  expect_identical(mult$z_score, mvn$z_score)
  expect_equal(mult$fwer_critical_value, mvn$fwer_critical_value, tolerance = 0.03)
})

test_that("hypo_test() is reproducible and leaves the caller's RNG alone", {
  toy <- toy_result()
  set.seed(42)
  before <- .Random.seed
  h1 <- hypo_test(toy, n_sim = 2000, seed = 7)
  expect_identical(.Random.seed, before)
  expect_identical(h1, hypo_test(toy, n_sim = 2000, seed = 7))
})

test_that("hypo_test() rejects invalid inputs", {
  toy <- toy_result()
  bad_len <- toy; bad_len$eif_mat <- toy$eif_mat[1:2]
  bad_est <- toy; bad_est$adj_association[1] <- NA

  expect_error(hypo_test(toy, n_sim = 0), "n_sim")
  expect_error(hypo_test(toy, fwer = 1.2), "fwer")
  expect_error(hypo_test(toy, seed = -1), "seed")
  expect_error(hypo_test(toy, method = "bonferroni"))
  expect_error(hypo_test(bad_len), "one entry per FC")
  expect_error(hypo_test(bad_est), "finite numeric vector")
})
