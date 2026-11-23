# Shared fast fit for the tests: one outcome (Frontoparietal) plus the
# seed-to-self column, a two-learner Super Learner library and log-normal GLM
# motion densities. Runs in about two seconds.
#
# The nuisance GLMs warn "fitted probabilities numerically 0 or 1 occurred" on
# the example data, so the warnings are suppressed here.
fast_moco <- function(..., Y = moco_example$Y[, 6:7]) {
  d <- moco_example
  args <- utils::modifyList(
    list(X = d$X, Z = d$Z, A = d$A, M = d$M, Y = Y,
         Delta_M = d$Delta_M, Delta_Y = d$Delta_Y,
         SL_library = c("SL.mean", "SL.glm"),
         density_method = "glm", cv_folds = 5, seed_rng = 1,
         hypo_n_sim = 2000),
    list(...),
    keep.null = TRUE
  )
  suppressWarnings(do.call(moco, args))
}
