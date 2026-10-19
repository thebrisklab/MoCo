toy_result <- function(n = 200, p = 3, seed = 1) {
  set.seed(seed)
  eif <- lapply(seq_len(p), function(j) cbind(eif_A0 = rnorm(n), eif_A1 = rnorm(n)))
  list(
    adj_association = stats::setNames(seq(0.3, 0, length.out = p), paste0("Y", seq_len(p))),
    eif_mat = eif,
    cov_mat = lapply(eif, function(x) stats::cov(x) / n)
  )
}

# Model fits on the small example data emit expected numerical warnings from
# the nuisance learners (e.g. rank-deficient or separated GLM fits).
quiet_moco <- function(...) suppressWarnings(moco(...))
