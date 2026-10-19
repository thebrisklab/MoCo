#' Motion-controlled estimation of group differences
#'
#' @description
#' `moco()` is the single high-level interface for MoCo. It estimates
#' group-specific means and their contrast while accounting for motion-related
#' selection. Conditional motion densities may be fitted with highly adaptive
#' lasso (HAL), log-normal generalized linear models (GLM), or generalized
#' additive models for location, scale, and shape (GAMLSS).
#'
#' @param X A data frame or matrix of baseline covariates, or `NULL` for no
#'   baseline covariates.
#' @param Z A data frame or matrix of post-exposure covariates used by the
#'   selection-bias adjustment, or `NULL`.
#' @param A A binary exposure or group indicator of length `n`.
#' @param M A numeric motion measure of length `n`. It must be observed
#'   whenever `Delta_Y = 1` or `Delta_M = 1`, and strictly positive for GLM and
#'   GAMLSS motion densities. Use `M = NULL` when motion is unavailable; the
#'   estimator then treats participants with `Delta_M = 0` as missing (as in
#'   Nebel et al., 2022) and no motion density is estimated.
#' @param Y A numeric outcome vector or an `n` by `p` numeric matrix. Rows with
#'   `Delta_Y = 0` may contain missing values. An entirely missing column is
#'   treated as a structural seed/self-correlation position and restored as
#'   `NA` in the output.
#' @param Delta_M A binary motion-inclusion vector. Supply `Delta_M` or
#'   `thresh`.
#' @param thresh Optional threshold used to define `Delta_M` as
#'   `M < thresh` (missing `M` gives `Delta_M = 0`). If supplied, it overrides
#'   `Delta_M`.
#' @param Delta_Y A binary outcome-observation and quality indicator.
#' @param SL_library Super Learner library used for nuisance regressions.
#'   Learners whose packages are not installed (for example `SL.xgboost`
#'   without \pkg{xgboost}) are dropped with a message.
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
#' @param gamlss_c_crit Convergence criterion for the GAMLSS global deviance,
#'   passed to [gamlss::gamlss.control()]. The default, `0.001`, is the
#'   `gamlss` default.
#' @param HAL_options Options passed to `haldensify` for HAL densities.
#' @param density_floor Strictly positive absolute lower bound used for motion
#'   densities appearing in denominators. The default,
#'   `sqrt(.Machine$double.eps)`, is a numerical safeguard against division by
#'   zero and matches the original MoCo implementation.
#' @param density_ratio_cap Optional positive upper bound for finite motion
#'   density ratios. The default, `Inf`, applies no cap and therefore leaves the
#'   paper estimand unchanged. A finite value truncates the density-ratio
#'   weights and therefore changes the estimand; use it only for sensitivity
#'   analyses of unstable density ratios.
#' @param support_rule Rule for where the target motion density
#'   `p(m | Delta_M = 1, A = 0, X)` is zero when it is evaluated at each
#'   participant's observed motion. `"delta_m"` (default) is the original MoCo
#'   rule: GLM and GAMLSS set it to zero for every participant with
#'   `Delta_M = 0`, and HAL uses the `haldensify` prediction, which is zero
#'   outside the training range. It is appropriate when `Delta_M` is a cutoff on
#'   `M` itself (for example `Delta_M = 1` if mean FD is below `thresh`).
#'   `"observed_range"` sets the density to zero only above the largest `M`
#'   among `Delta_M = 1` participants in the full sample; use it when `Delta_M`
#'   is defined by scrubbing or a retained-frame count, so that participants
#'   with `Delta_M = 0` can have tolerable motion.
#' @param cross_fit Logical indicating whether to use cross-fitting.
#' @param cv_folds Number of cross-fitting folds.
#' @param test Logical indicating whether to run simultaneous EIF tests.
#' @param fwer Numeric vector of family-wise error rates.
#' @param seed_rgn Integer vector defining repeated nuisance-model and
#'   cross-fitting splits, not optimizer initializations. Seeds are used
#'   internally and the caller's random-number state is restored on exit.
#' @param test_seed Seed used by [hypo_test()].
#' @param test_n_sim Number of multivariate-normal draws used by [hypo_test()].
#' @param test_chunk_size Maximum number of simultaneous-test draws generated
#'   at once. The default, `test_n_sim`, reproduces the original MoCo critical
#'   values; smaller values reduce peak memory use.
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
#' influence functions are averaged across seeds, and the covariance,
#' z-statistics, and simultaneous critical values are then recomputed from the
#' averaged influence functions. This replaces the averaging of z-statistics
#' and critical values across runs used in earlier versions of MoCo; the point
#' estimates are the same under both approaches, and a single seed gives
#' identical results.
#'
#' @return A list with components:
#' \describe{
#'   \item{est}{Data frame with two rows (`est_A0`, `est_A1`) and one column
#'     per outcome: motion-controlled group means.}
#'   \item{adj_association}{Named numeric vector of `A = 1` minus `A = 0`
#'     contrasts.}
#'   \item{density_method}{Character vector naming the `pMX` and `pMXZ`
#'     density methods.}
#'   \item{density_controls}{List recording the absolute denominator floor and
#'     optional density-ratio cap used in the fit.}
#'   \item{gamlss_selection_runs}{(GAMLSS only) List of fitted GAMLSS
#'     specifications for each seed.}
#'   \item{z_score}{(`test = TRUE`) Named vector of EIF-based z-statistics.}
#'   \item{conf_band}{(`test = TRUE`) Simultaneous critical value for each
#'     requested FWER.}
#'   \item{significant_regions}{(`test = TRUE`) Logical vector (one FWER) or
#'     matrix (one row per FWER) of simultaneous-test decisions.}
#' }
#' Structurally missing outcome columns are returned as `NA`.
#' @seealso [hypo_test()]
#' @examples
#' data(data)
#'
#' \donttest{
#' # A fast illustration: log-normal GLM motion densities, a small learner
#' # library, two outcomes and two cross-fitting folds.
#' fit <- moco(
#'   X = data$X, Z = data$Z, A = data$A, M = data$M,
#'   Y = data$Y[, 1:2],
#'   Delta_M = data$Delta_M, Delta_Y = data$Delta_Y,
#'   SL_library = c("SL.mean", "SL.glm"),
#'   pMX_method = "GLM",
#'   glm_formula = list(pMX = ".", pMXZ = "."),
#'   cv_folds = 2,
#'   test_n_sim = 2000
#' )
#' fit$est
#' fit$adj_association
#' fit$significant_regions
#' }
#' @import SuperLearner
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
    gamlss_c_crit = 0.001,
    HAL_options = list(
      max_degree = 3,
      lambda_seq = exp(seq(-1, -10, length = 100)),
      num_knots = c(1000, 500, 250)
    ),
    density_floor = sqrt(.Machine$double.eps),
    density_ratio_cap = Inf,
    support_rule = c("delta_m", "observed_range"),
    cross_fit = TRUE,
    cv_folds = 5,
    test = TRUE,
    fwer = 0.05,
    seed_rgn = 1,
    test_seed = 1,
    test_n_sim = 100000L,
    test_chunk_size = test_n_sim,
    ...
) {
  # Seeds are set internally for reproducible nuisance fits; restore the
  # caller's random-number state on exit.
  withr::local_preserve_seed()

  SL_library <- .moco_available_learners(SL_library)
  SL_library_customize <- lapply(
    stats::setNames(names(SL_library_customize), names(SL_library_customize)),
    function(component) {
      .moco_available_learners(
        SL_library_customize[[component]],
        paste0("SL_library_customize$", component)
      )
    }
  )

  support_rule <- match.arg(support_rule)
  M_available <- !is.null(M)

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
  if (!M_available && !is.null(thresh)) {
    stop("thresh requires M; supply Delta_M when M = NULL.")
  }
  if (!is.null(Delta_M) && !is.null(thresh)) {
    warning("Both Delta_M and thresh were supplied; Delta_M is recomputed as M < thresh.")
  }
  density_controls <- .moco_validate_density_controls(
    density_floor, density_ratio_cap
  )

  if (!is.null(X)) X <- as.data.frame(X)
  if (!is.null(Z)) Z <- as.data.frame(Z)
  if (anyNA(A) || anyNA(X) || anyNA(Z)) {
    stop("A, X, and Z must not contain missing values.")
  }

  Y <- as.matrix(Y)
  if (!is.numeric(Y)) {
    stop("Y must be numeric.")
  }
  n <- nrow(Y)
  if (
    length(A) != n || (M_available && length(M) != n) ||
      (!is.null(X) && nrow(X) != n) || (!is.null(Z) && nrow(Z) != n) ||
      length(Delta_Y) != n
  ) {
    stop("A, M, X, Z, Delta_Y, and Y must have the same number of observations.")
  }
  if (anyNA(Delta_Y) || any(!Delta_Y %in% c(0, 1))) {
    stop("Delta_Y must contain only 0 and 1 with no missing values.")
  }
  if (M_available && !is.numeric(M)) {
    stop("M must be numeric.")
  }
  if (M_available && any(Delta_Y == 1 & is.na(M))) {
    stop("M must be observed for every row with Delta_Y = 1.")
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

  if (M_available && !is.null(Delta_M) && is.null(thresh) &&
      any(Delta_M == 1 & is.na(M))) {
    stop("M must be observed for every row with Delta_M = 1.")
  }

  density_methods <- if (M_available) {
    c(pMX = pMX_method, pMXZ = pMXZ_method)
  } else {
    c(pMX = "none", pMXZ = "none")
  }
  if (any(density_methods %in% c("GLM", "GAMLSS")) &&
      any(!is.na(M) & M <= 0)) {
    stop("GLM and GAMLSS motion densities require strictly positive M.")
  }
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
      density_floor = density_controls$density_floor,
      density_ratio_cap = density_controls$density_ratio_cap,
      support_rule = support_rule,
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
          n.cyc = n.cyc,
          gamlss_c_crit = gamlss_c_crit
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
    density_method = density_methods,
    density_controls = density_controls,
    support_rule = support_rule
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
