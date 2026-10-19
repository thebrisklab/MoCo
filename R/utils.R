#' Evaluate EIF for stochastic interventional (in)direct effects
#'
#' @param a Either 0 or 1.
#' @param A A binary exposure or group-indicator vector of length `n`.
#' @param Delta_Y Binary indicator that the outcome is observed and usable.
#' @param Delta_M Binary indicator that motion meets the inclusion criterion.
#' @param Y A vector of continuous outcome of interest.
#' @param gA Estimate of P(A = 1 | X_i), i = 1,...,n.
#' @param gDM Estimate of P(Delta_M = 1 | A_i, X_i), i = 1,...,n.
#' @param gDY_AXZ Estimate of P(Delta_Y = 1 | A_i, X_i, Z_i).
#' @param eta_AXZ Estimate of E(mu_pseudo_A | A_i, Delta_i, X_i).
#' @param eta_AXM Estimate of E(mu_pseudo_A_star | A_i, M_i, X_i).
#' @param xi_AX Estimate of E(eta_AXZ | A_i, Delta_i, X_i).
#' @param mu Estimate of E(Y| A_i, M_i, X_i, Z_i).
#' @param pMXD Estimate of p(M | A = 0, Delta_M = 1, X_i).
#' @param pMXZ Estimate of p(M | A_i, X_i, Z_i).
#'
#' @return An n-length vector of the estimated EIF evaluated on the observations.
#' @noRd

make_full_data_eif <- function(a, A, Delta_Y, Delta_M, Y, gA, gDM, gDY_AXZ,
                               eta_AXZ, eta_AXM, xi_AX, mu, pMXD, pMXZ,
                               density_ratio_cap = Inf){
  n <- length(A)
  # Each term is non-zero only for the participants who contribute to it;
  # elsewhere its inputs may legitimately be missing.
  use_p1 <- Delta_Y == 1 & A == a
  use_p2 <- A == a
  use_p3 <- A == 0 & Delta_M == 1

  p1 <- p2 <- p3 <- numeric(n)
  density_ratio <- .moco_cap_density_ratio(pMXD / pMXZ, density_ratio_cap)
  p1[use_p1] <- (density_ratio * (Y - mu) / (gDY_AXZ * gA))[use_p1]
  p2[use_p2] <- ((eta_AXZ - xi_AX) / gA)[use_p2]
  gA_0 <- (1 - gA) * (a == 1) + gA * (a == 0)
  p3[use_p3] <- ((eta_AXM - xi_AX) / (gA_0 * gDM))[use_p3]
  p4 <- xi_AX - mean(xi_AX)

  .moco_check_eif_terms(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
  p1 + p2 + p3 + p4
}

make_full_data_eif_easy <- function(a, A, Delta_Y, Y, gA, gDY_AXZ, xi_AX, mu) {
  n <- length(A)
  use_p1 <- Delta_Y == 1 & A == a
  use_p2 <- A == a

  p1 <- p2 <- numeric(n)
  p1[use_p1] <- ((Y - mu) / (gDY_AXZ * gA))[use_p1]
  p2[use_p2] <- ((mu - xi_AX) / gA)[use_p2]
  p3 <- xi_AX - mean(xi_AX)

  .moco_check_eif_terms(list(p1 = p1, p2 = p2, p3 = p3))
  p1 + p2 + p3
}

# Stop with an informative error instead of silently treating missing or
# infinite nuisance predictions as zero.
.moco_check_eif_terms <- function(terms) {
  bad <- names(terms)[!vapply(terms, function(x) all(is.finite(x)), logical(1))]
  if (length(bad) > 0L) {
    stop(
      "Non-finite values in EIF term(s) ", paste(bad, collapse = ", "),
      ". Check the nuisance fits (missing motion values, zero densities or ",
      "zero estimated probabilities).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Average repeated one-step fits at the estimator and participant-level EIF
# stages, then recompute covariance. Averaging z statistics or critical values
# directly does not represent the variance of the averaged estimator.
.moco_combine_one_step <- function(results) {
  if (!is.list(results) || length(results) == 0L) {
    stop("results must be a non-empty list of one-step fit objects.")
  }

  required <- c("est", "eif_mat")
  for (i in seq_along(results)) {
    missing_components <- setdiff(required, names(results[[i]]))
    if (length(missing_components) > 0L) {
      stop(
        "Seed fit ", i, " is missing: ",
        paste(missing_components, collapse = ", "), "."
      )
    }
  }

  reference_est <- as.matrix(results[[1L]]$est)
  if (nrow(reference_est) != 2L || ncol(reference_est) < 1L ||
      any(!is.finite(reference_est))) {
    stop("Each seed fit must contain a finite 2 by p est matrix.")
  }
  n_features <- ncol(reference_est)
  if (!is.list(results[[1L]]$eif_mat) ||
      length(results[[1L]]$eif_mat) != n_features) {
    stop("Each seed fit must contain one EIF matrix per outcome.")
  }
  reference_eif <- as.matrix(results[[1L]]$eif_mat[[1L]])
  if (ncol(reference_eif) != 2L || nrow(reference_eif) < 2L) {
    stop("Each EIF must be an n by 2 matrix.")
  }
  n_obs <- nrow(reference_eif)

  normalized <- lapply(seq_along(results), function(i) {
    est_i <- as.matrix(results[[i]]$est)
    if (!identical(dim(est_i), c(2L, n_features)) ||
        any(!is.finite(est_i))) {
      stop("All seed fits must have finite est matrices with identical dimensions.")
    }
    if (!is.list(results[[i]]$eif_mat) ||
        length(results[[i]]$eif_mat) != n_features) {
      stop("All seed fits must contain the same number of EIF matrices.")
    }
    eif_i <- lapply(results[[i]]$eif_mat, function(x) {
      x <- as.matrix(x)
      if (!identical(dim(x), c(n_obs, 2L)) || any(!is.finite(x))) {
        stop(
          "All seed EIF matrices must be finite, have dimension n by 2, ",
          "and retain the same participant order."
        )
      }
      x
    })
    list(est = est_i, eif_mat = eif_i)
  })

  est <- Reduce(`+`, lapply(normalized, `[[`, "est")) / length(normalized)
  eif_mat <- lapply(seq_len(n_features), function(j) {
    Reduce(`+`, lapply(normalized, function(x) x$eif_mat[[j]])) /
      length(normalized)
  })
  cov_mat <- lapply(eif_mat, function(x) stats::cov(x) / n_obs)

  list(
    est = est,
    adj_association = est[2L, ] - est[1L, ],
    eif_mat = eif_mat,
    cov_mat = cov_mat
  )
}

#' Temporary fix for convex combination method mean squared error
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another
#' @noRd
tmp_method.CC_LS <- function() {
  computeCoef <- function(Z, Y, libraryNames, verbose, obsWeights,
                          errorsInLibrary = NULL, ...) {
    cvRisk <- apply(Z, 2, function(x) {
      mean(obsWeights * (x -
                           Y)^2)
    })
    names(cvRisk) <- libraryNames
    compute <- function(x, y, wt = rep(1, length(y))) {
      wX <- sqrt(wt) * x
      wY <- sqrt(wt) * y
      D <- crossprod(wX)
      d <- crossprod(wX, wY)
      A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
      bvec <- c(1, rep(0, ncol(wX)))
      fit <- tryCatch(
        {
          quadprog::solve.QP(
            Dmat = D, dvec = d, Amat = A,
            bvec = bvec, meq = 1
          )
        },
        error = function(e) {
          out <- list()
          class(out) <- "error"
          out
        }
      )
      invisible(fit)
    }
    modZ <- Z
    naCols <- which(apply(Z, 2, function(z) {
      all(z == 0)
    }))
    anyNACols <- length(naCols) > 0
    if (anyNACols) {
      warning(paste0(
        paste0(libraryNames[naCols], collapse = ", "),
        " have NAs.", "Removing from super learner."
      ))
    }
    tol <- 4
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    # if (anyDupCols) {
    #   warning(paste0(
    #     paste0(libraryNames[dupCols], collapse = ", "),
    #     " are duplicates of previous learners.", " Removing from super learner."
    #   ))
    # }
    if (anyDupCols | anyNACols) {
      rmCols <- unique(c(naCols, dupCols))
      modZ <- Z[, -rmCols, drop = FALSE]
    }
    fit <- compute(x = modZ, y = Y, wt = obsWeights)
    if (!inherits(fit, "error")) {
      coef <- fit$solution
    } else {
      coef <- rep(0, ncol(Z))
      coef[which.min(cvRisk)] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (!inherits(fit, "error")) {
      if (anyDupCols | anyNACols) {
        full_coef <- rep(0, ncol(Z))
        full_coef[-rmCols] <- coef
        coef <- full_coef
      }
      coef[coef < 1e-04] <- 0
      coef <- coef / sum(coef)
    }
    if (!sum(coef) > 0) {
      warning("All algorithms have zero weight", call. = FALSE)
    }
    list(cvRisk = cvRisk, coef = coef, optimizer = fit)
  }
  computePred <- function(predY, coef, ...) {
    predY %*% matrix(coef)
  }
  out <- list(
    require = "quadprog", computeCoef = computeCoef,
    computePred = computePred
  )
  invisible(out)
}


#' Temporary fix for convex combination method negative log-likelihood loss
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another.
#' Note that because of the way \code{SuperLearner} is structure, one needs to
#' install the optimization software separately.
#' @noRd
tmp_method.CC_nloglik <- function() {
  computePred <- function(predY, coef, control, ...) {
    if (sum(coef != 0) == 0) {
      stop("All metalearner coefficients are zero, cannot compute prediction.")
    }
    stats::plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*%
                    matrix(coef[coef != 0]))
  }
  computeCoef <- function(Z, Y, libraryNames, obsWeights, control,
                          verbose, ...) {
    tol <- 4
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    modZ <- Z
    if (anyDupCols) {
      # warning(paste0(
      #   paste0(libraryNames[dupCols], collapse = ", "),
      #   " are duplicates of previous learners.", " Removing from super learner."
      # ))
      modZ <- modZ[, -dupCols, drop = FALSE]
    }
    modlogitZ <- trimLogit(modZ, control$trimLogit)
    logitZ <- trimLogit(Z, control$trimLogit)
    cvRisk <- apply(logitZ, 2, function(x) {
      -sum(2 * obsWeights *
             ifelse(Y, stats::plogis(x, log.p = TRUE), stats::plogis(x,
                                                                     log.p = TRUE,
                                                                     lower.tail = FALSE
             )))
    })
    names(cvRisk) <- libraryNames
    obj_and_grad <- function(y, x, w = NULL) {
      y <- y
      x <- x
      function(beta) {
        xB <- x %*% cbind(beta)
        loglik <- y * stats::plogis(xB, log.p = TRUE) + (1 -
                                                           y) * stats::plogis(xB, log.p = TRUE, lower.tail = FALSE)
        if (!is.null(w)) {
          loglik <- loglik * w
        }
        obj <- -2 * sum(loglik)
        p <- stats::plogis(xB)
        grad <- if (is.null(w)) {
          2 * crossprod(x, cbind(p - y))
        } else {
          2 * crossprod(x, w * cbind(p - y))
        }
        list(objective = obj, gradient = grad)
      }
    }
    lower_bounds <- rep(0, ncol(modZ))
    upper_bounds <- rep(1, ncol(modZ))
    if (anyNA(cvRisk)) {
      upper_bounds[is.na(cvRisk)] <- 0
    }
    r <- tryCatch(
      {
        nloptr::nloptr(
          x0 = rep(1 / ncol(modZ), ncol(modZ)),
          eval_f = obj_and_grad(Y, modlogitZ, w = obsWeights), 
          lb = lower_bounds,
          ub = upper_bounds, eval_g_eq = function(beta) {
            (sum(beta) -
               1)
          }, eval_jac_g_eq = function(beta) rep(1, length(beta)),
          opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08)
        )
      },
      error = function(e) {
        out <- list()
        class(out) <- "error"
        out
      }
    )
    optimizer_failed <- inherits(r, "error")
    if (!optimizer_failed && (r$status < 1 || r$status > 4)) {
      warning(r$message)
    }
    if (!optimizer_failed) {
      coef <- r$solution
    } else {
      coef <- rep(0, ncol(Z))
      coef[which.min(cvRisk)] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    # The fallback coefficients already have one entry per learner, so only
    # re-insert the removed duplicate columns after a successful optimization.
    if (!optimizer_failed && anyDupCols) {
      full_coef <- rep(0, ncol(Z))
      full_coef[-dupCols] <- coef
      coef <- full_coef
    }
    coef[coef < 1e-04] <- 0
    coef <- coef / sum(coef)
    out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
    return(out)
  }
  list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}

# Packages needed by Super Learner wrappers whose packages are listed in
# Suggests rather than Imports.
.moco_learner_packages <- c(
  SL.earth = "earth",
  SL.gam = "gam",
  SL.glmnet = "glmnet",
  SL.ranger = "ranger",
  SL.xgboost = "xgboost"
)

# Drop learners whose package is not installed so that the default library
# works on systems without the suggested packages.
.moco_available_learners <- function(library, argument = "SL_library") {
  if (is.null(library) || !is.character(library)) return(library)
  needed <- .moco_learner_packages[intersect(library, names(.moco_learner_packages))]
  missing_pkg <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkg) > 0L) {
    message(
      "Dropping learners from ", argument, " because their packages are not installed: ",
      paste0(names(missing_pkg), " (", missing_pkg, ")", collapse = ", "), "."
    )
    library <- setdiff(library, names(missing_pkg))
  }
  if (length(library) == 0L) {
    stop(argument, " has no learners left after removing unavailable ones.")
  }
  library
}

# Largest observed motion value among participants who meet the motion
# inclusion criterion (NA if there are none).
.moco_support_upper <- function(M, Delta_M) {
  included <- M[!is.na(Delta_M) & Delta_M == 1 & is.finite(M)]
  if (length(included) == 0L) return(NA_real_)
  max(included)
}

# Apply the support rule to the target density p(m | Delta_M = 1, A = 0, X)
# evaluated at each participant's observed M.
# - "delta_m" (original MoCo): zero for every Delta_M = 0 participant when the
#   density is GLM or GAMLSS; HAL predictions are used as returned by
#   haldensify, which are zero outside the training range. Appropriate when
#   Delta_M is a cutoff on M itself (Delta_M = 1 if M < thresh).
# - "observed_range": zero only above the largest observed M among
#   Delta_M = 1 participants (full sample when cross-fitting). Appropriate when
#   Delta_M is defined by scrubbing or a retained-frame count.
.moco_apply_support_rule <- function(density, M, Delta_M, support_rule,
                                     support_upper = NULL, hal = FALSE) {
  if (identical(support_rule, "delta_m")) {
    if (!isTRUE(hal)) density[Delta_M == 0] <- 0
    return(density)
  }
  if (!identical(support_rule, "observed_range")) {
    stop("support_rule must be \"delta_m\" or \"observed_range\".")
  }
  if (is.null(support_upper)) support_upper <- .moco_support_upper(M, Delta_M)
  if (length(support_upper) != 1L || !is.finite(support_upper)) {
    stop("No finite motion values among observations with Delta_M = 1.")
  }
  density[!is.na(M) & M > support_upper] <- 0
  density
}

# Bound a density used as a denominator away from zero using an absolute floor.
# This is a numerical guard only; it is not relative to the empirical density
# distribution and therefore does not intentionally trim the target estimand.
.moco_floor_density <- function(density, density_floor = sqrt(.Machine$double.eps)) {
  density[!is.na(density) & density < density_floor] <- density_floor
  density
}

.moco_cap_density_ratio <- function(ratio, density_ratio_cap = Inf) {
  if (is.finite(density_ratio_cap)) {
    ratio[!is.na(ratio) & ratio > density_ratio_cap] <- density_ratio_cap
  }
  ratio
}

.moco_validate_density_controls <- function(
    density_floor = sqrt(.Machine$double.eps),
    density_ratio_cap = Inf
) {
  if (!is.numeric(density_floor) || length(density_floor) != 1L ||
      is.na(density_floor) || !is.finite(density_floor) ||
      density_floor <= 0) {
    stop("density_floor must be one strictly positive finite number.")
  }
  if (!is.numeric(density_ratio_cap) || length(density_ratio_cap) != 1L ||
      is.na(density_ratio_cap) || density_ratio_cap <= 0) {
    stop("density_ratio_cap must be one positive number or Inf.")
  }
  list(
    density_floor = as.numeric(density_floor),
    density_ratio_cap = as.numeric(density_ratio_cap)
  )
}
