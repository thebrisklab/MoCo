#' Motion-controlled estimation of group differences
#'
#' @description
#' `moco()` is the single high-level interface for MoCo. It estimates
#' group-specific means and their contrast while accounting for motion-related
#' selection. Conditional motion densities may be fitted with highly adaptive
#' lasso (HAL), log-normal generalized linear models (GLM), or generalized
#' additive models for location, scale, and shape (GAMLSS).
#'
#' @param X A data frame or matrix of baseline covariates.
#' @param Z A data frame or matrix of post-exposure covariates used by the
#'   selection-bias adjustment.
#' @param A A binary exposure or group indicator of length `n`.
#' @param M A numeric motion measure of length `n`.
#' @param Y A numeric outcome vector or an `n` by `p` numeric matrix. Rows with
#'   `Delta_Y = 0` may contain missing values. An entirely missing column is
#'   treated as a structural seed/self-correlation position and restored as
#'   `NA` in the output.
#' @param Delta_M A binary motion-inclusion vector. Supply `Delta_M` or
#'   `thresh`.
#' @param thresh Optional threshold used to define `Delta_M` from `M`.
#' @param Delta_Y A binary outcome-observation and quality indicator.
#' @param SL_library Super Learner library used for nuisance regressions.
#' @param SL_library_customize Named list optionally specifying separate
#'   libraries for `gA`, `gDM`, `gDY_AX`, `gDY_AXZ`, `mu_AMXZ`, `eta_AXZ`,
#'   `eta_AXM`, and `xi_AX`.
#' @param glm_formula Named list of right-hand-side formulas for GLM nuisance
#'   regressions. For a GLM motion density, provide `pMX` and/or `pMXZ`, for
#'   example `list(pMX = ".", pMXZ = ".")`.
#' @param pMX_method Method for densities conditional on `A` and `X`: `"HAL"`,
#'   `"GLM"`, or `"GAMLSS"`.
#' @param pMXZ_method Method for densities conditional on `A`, `X`, and `Z`.
#'   If `NULL`, it uses `pMX_method`.
#' @param gamlss_formula Reserved for compatibility. Current GAMLSS formulas
#'   are built from `gamlss_continuous_X` and `gamlss_continuous_Z`.
#' @param gamlss_family GAMLSS family. Currently only generalized gamma
#'   (`"GG"`) is supported.
#' @param gamlss_optimizer GAMLSS algorithm: `"RS"`, `"CG"`, or `"mixed"`.
#' @param GAMLSS_BIC_select Reserved for compatibility; model-structure
#'   selection is not performed.
#' @param gamlss_continuous_X Continuous columns of `X` that enter GAMLSS
#'   location and scale models through [gamlss::pb()].
#' @param gamlss_continuous_Z Continuous columns of `Z` that enter GAMLSS
#'   location and scale models through [gamlss::pb()].
#' @param gamlss_bic_candidates Reserved for compatibility.
#' @param gamlss_bic_trace Logical controlling GAMLSS fitting output.
#' @param n.cyc Maximum number of GAMLSS fitting cycles.
#' @param HAL_options Options passed to `haldensify` for HAL densities.
#' @param cross_fit Logical indicating whether to use cross-fitting.
#' @param cv_folds Number of cross-fitting folds.
#' @param test Logical indicating whether to run simultaneous EIF tests.
#' @param fwer Numeric vector of family-wise error rates.
#' @param seed_rgn Integer vector defining repeated nuisance-model and
#'   cross-fitting splits, not optimizer initializations.
#' @param test_seed Seed used by [hypo_test()].
#' @param test_n_sim Number of multivariate-normal draws used by [hypo_test()].
#' @param test_chunk_size Maximum number of simultaneous-test draws generated
#'   at once. Smaller values reduce peak memory use.
#' @param ... Additional named arguments passed to the selected fitting engine.
#'
#' @details
#' `pMX_method` and `pMXZ_method` may be selected independently, although most
#' analyses use the same method for both. GLM densities are log-normal. GAMLSS
#' densities are generalized gamma, with continuous covariates entering the
#' location and scale models through penalized splines, other covariates
#' entering linearly, and a constant shape parameter.
#'
#' For repeated `seed_rgn` values, estimates and participant-level efficient
#' influence functions are averaged before covariance, z-statistics, and
#' simultaneous critical values are recomputed.
#'
#' @return A list with `est`, `adj_association`, and `density_method`. GAMLSS
#'   fits also return `gamlss_selection_runs`. When `test = TRUE`, the list
#'   additionally contains `z_score`, `conf_band`, and `significant_regions`.
#' @seealso [hypo_test()]
#' @import SuperLearner
#' @import haldensify
#' @import MASS
#' @export
moco <- function(
    X,
    Z,
    A,
    M,
    Y,
    Delta_M = NULL,
    thresh = NULL,
    Delta_Y,
    SL_library = c(
      "SL.earth", "SL.glmnet", "SL.gam", "SL.glm",
      "SL.glm.interaction", "SL.step", "SL.step.interaction",
      "SL.xgboost", "SL.ranger", "SL.mean"
    ),
    SL_library_customize = list(
      gA = NULL, gDM = NULL, gDY_AX = NULL, gDY_AXZ = NULL,
      mu_AMXZ = NULL, eta_AXZ = NULL, eta_AXM = NULL, xi_AX = NULL
    ),
    glm_formula = list(
      gA = NULL, gDM = NULL, gDY_AX = NULL, gDY_AXZ = NULL,
      mu_AMXZ = NULL, eta_AXZ = NULL, eta_AXM = NULL, xi_AX = NULL,
      pMX = NULL, pMXZ = NULL
    ),
    pMX_method = c("HAL", "GLM", "GAMLSS"),
    pMXZ_method = NULL,
    gamlss_formula = list(
      pMX_mu = NULL, pMX_sigma = ~ 1, pMX_nu = ~ 1,
      pMXZ_mu = NULL, pMXZ_sigma = ~ 1, pMXZ_nu = ~ 1
    ),
    gamlss_family = "GG",
    gamlss_optimizer = "RS",
    GAMLSS_BIC_select = FALSE,
    gamlss_continuous_X = character(0),
    gamlss_continuous_Z = character(0),
    gamlss_bic_candidates = c("linear", "pb_mu", "pb_mu_sigma"),
    gamlss_bic_trace = FALSE,
    n.cyc = 300,
    HAL_options = list(
      max_degree = 3,
      lambda_seq = exp(seq(-1, -10, length = 100)),
      num_knots = c(1000, 500, 250)
    ),
    cross_fit = TRUE,
    cv_folds = 5,
    test = TRUE,
    fwer = 0.05,
    seed_rgn = 1,
    test_seed = 1,
    test_n_sim = 100000L,
    test_chunk_size = 25000L,
    ...
) {
  pMX_method <- .moco_match_density_method(pMX_method, "pMX_method")
  if (is.null(pMXZ_method)) {
    pMXZ_method <- pMX_method
  } else {
    pMXZ_method <- .moco_match_density_method(pMXZ_method, "pMXZ_method")
  }

  if (length(seed_rgn) == 0L || anyNA(seed_rgn)) {
    stop("seed_rgn must contain at least one non-missing seed.")
  }
  cv_folds <- as.integer(cv_folds)
  if (length(cv_folds) != 1L || is.na(cv_folds) || cv_folds < 2L) {
    stop("cv_folds must be one integer greater than or equal to 2.")
  }
  if (is.null(Delta_M) && is.null(thresh)) {
    stop("Supply either Delta_M or thresh.")
  }

  X <- as.data.frame(X)
  Z <- as.data.frame(Z)
  if (anyNA(A) || anyNA(X) || anyNA(Z)) {
    stop("A, X, and Z must not contain missing values.")
  }

  Y <- as.matrix(Y)
  if (!is.numeric(Y)) {
    stop("Y must be numeric.")
  }
  n <- nrow(Y)
  if (
    length(A) != n || length(M) != n || nrow(X) != n || nrow(Z) != n ||
      length(Delta_Y) != n
  ) {
    stop("A, M, X, Z, Delta_Y, and Y must have the same number of observations.")
  }
  if (anyNA(Delta_Y) || any(!Delta_Y %in% c(0, 1))) {
    stop("Delta_Y must contain only 0 and 1 with no missing values.")
  }
  if (is.null(thresh)) {
    if (length(Delta_M) != n || anyNA(Delta_M) || any(!Delta_M %in% c(0, 1))) {
      stop("Delta_M must contain only 0 and 1 with no missing values.")
    }
  } else if (length(thresh) != 1L || is.na(thresh) || !is.finite(thresh)) {
    stop("thresh must be one finite numeric value.")
  }

  outcome_names <- colnames(Y)
  if (is.null(outcome_names)) outcome_names <- paste0("Y", seq_len(ncol(Y)))
  structural_positions <- which(colSums(!is.na(Y)) == 0L)
  analysis_positions <- setdiff(seq_len(ncol(Y)), structural_positions)
  if (length(analysis_positions) == 0L) {
    stop("Y must contain at least one outcome that is not entirely missing.")
  }
  Y_analysis <- Y[, analysis_positions, drop = FALSE]
  if (any(!is.finite(Y_analysis[Delta_Y == 1, , drop = FALSE]))) {
    stop("Y must be finite for every row with Delta_Y = 1.")
  }

  density_methods <- c(pMX = pMX_method, pMXZ = pMXZ_method)
  glm_components <- names(density_methods)[density_methods == "GLM"]
  missing_glm_formula <- glm_components[vapply(
    glm_components,
    function(component) is.null(glm_formula[[component]]),
    logical(1)
  )]
  if (length(missing_glm_formula) > 0L) {
    stop(
      "GLM density selected but glm_formula is missing: ",
      paste(missing_glm_formula, collapse = ", "), "."
    )
  }

  use_gamlss <- any(density_methods == "GAMLSS")
  if (use_gamlss) {
    if (!identical(gamlss_family, "GG")) {
      stop("The current GAMLSS density implementation supports gamlss_family = \"GG\" only.")
    }
    gamlss_optimizer <- .moco_match_gamlss_optimizer(gamlss_optimizer)
    missing_cx <- setdiff(gamlss_continuous_X, names(X))
    missing_cz <- setdiff(gamlss_continuous_Z, names(Z))
    if (length(missing_cx) > 0L) {
      stop("gamlss_continuous_X not found in X: ", paste(missing_cx, collapse = ", "))
    }
    if (length(missing_cz) > 0L) {
      stop("gamlss_continuous_Z not found in Z: ", paste(missing_cz, collapse = ", "))
    }
    if (isTRUE(GAMLSS_BIC_select)) {
      warning("GAMLSS_BIC_select is ignored; using the default pb() specification.")
    }
  }

  HAL_pMX <- identical(pMX_method, "HAL")
  HAL_pMXZ <- identical(pMXZ_method, "HAL")
  GAMLSS_pMX <- identical(pMX_method, "GAMLSS")
  GAMLSS_pMXZ <- identical(pMXZ_method, "GAMLSS")
  extra_args <- list(...)
  if (length(extra_args) > 0L &&
      (is.null(names(extra_args)) || any(!nzchar(names(extra_args))))) {
    stop("All arguments supplied through ... must be named.")
  }

  fits <- vector("list", length(seed_rgn))
  gamlss_selection_runs <- vector("list", length(seed_rgn))

  for (i in seq_along(seed_rgn)) {
    common_args <- list(
      X = X, Z = Z, A = A, M = M, Y = Y_analysis,
      Delta_M = Delta_M, thresh = thresh, Delta_Y = Delta_Y,
      SL_library = SL_library,
      SL_library_customize = SL_library_customize,
      glm_formula = glm_formula,
      HAL_pMX = HAL_pMX, HAL_pMXZ = HAL_pMXZ,
      HAL_options = HAL_options,
      seed = seed_rgn[i]
    )

    if (use_gamlss) {
      engine <- if (isTRUE(cross_fit)) one_step_cross_gamlss else one_step_gamlss
      engine_args <- c(
        common_args,
        list(
          GAMLSS_pMX = GAMLSS_pMX,
          GAMLSS_pMXZ = GAMLSS_pMXZ,
          gamlss_formula = gamlss_formula,
          gamlss_family = gamlss_family,
          gamlss_optimizer = gamlss_optimizer,
          GAMLSS_BIC_select = GAMLSS_BIC_select,
          gamlss_continuous_X = gamlss_continuous_X,
          gamlss_continuous_Z = gamlss_continuous_Z,
          gamlss_bic_candidates = gamlss_bic_candidates,
          gamlss_bic_trace = gamlss_bic_trace,
          n.cyc = n.cyc
        )
      )
    } else {
      engine <- if (isTRUE(cross_fit)) one_step_cross else one_step
      engine_args <- common_args
    }
    if (isTRUE(cross_fit)) engine_args$cv_folds <- cv_folds
    duplicate_dots <- intersect(names(engine_args), names(extra_args))
    if (length(duplicate_dots) > 0L) {
      stop("Arguments supplied more than once: ", paste(duplicate_dots, collapse = ", "), ".")
    }

    result <- do.call(engine, c(engine_args, extra_args))
    if (use_gamlss) {
      if (!is.null(result$gamlss_selection_by_fold)) {
        gamlss_selection_runs[[i]] <- result$gamlss_selection_by_fold
      } else if (!is.null(result$gamlss_selection)) {
        gamlss_selection_runs[[i]] <- result$gamlss_selection
      }
    }
    fits[[i]] <- result
  }

  combined <- .moco_combine_one_step(fits)
  n_outcomes <- ncol(Y)
  est_full <- matrix(
    NA_real_, 2L, n_outcomes,
    dimnames = list(c("est_A0", "est_A1"), outcome_names)
  )
  est_full[, analysis_positions] <- combined$est
  adj_full <- stats::setNames(rep(NA_real_, n_outcomes), outcome_names)
  adj_full[analysis_positions] <- combined$adj_association

  output <- list(
    est = as.data.frame(est_full),
    adj_association = adj_full,
    density_method = density_methods
  )
  if (use_gamlss) output$gamlss_selection_runs <- gamlss_selection_runs

  if (isTRUE(test)) {
    hypo <- hypo_test(
      combined,
      fwer = fwer,
      seed = test_seed,
      n_sim = test_n_sim,
      chunk_size = test_chunk_size
    )
    z_full <- stats::setNames(rep(NA_real_, n_outcomes), outcome_names)
    z_full[analysis_positions] <- hypo$z_score
    significant_full <- matrix(
      NA,
      nrow = nrow(hypo$significant_regions),
      ncol = n_outcomes,
      dimnames = list(rownames(hypo$significant_regions), outcome_names)
    )
    significant_full[, analysis_positions] <- hypo$significant_regions
    if (length(fwer) == 1L) {
      significant_full <- stats::setNames(as.logical(significant_full[1L, ]), outcome_names)
    }
    output <- c(
      output,
      list(
        z_score = z_full,
        conf_band = hypo$conf_band,
        significant_regions = significant_full
      )
    )
  }

  output
}

.moco_match_density_method <- function(method, argument) {
  choices <- c("HAL", "GLM", "GAMLSS")
  if (length(method) > 1L) method <- method[1L]
  if (length(method) != 1L || is.na(method)) {
    stop(argument, " must be one of: ", paste(choices, collapse = ", "), ".")
  }
  method <- toupper(as.character(method))
  if (!method %in% choices) {
    stop(argument, " must be one of: ", paste(choices, collapse = ", "), ".")
  }
  method
}
