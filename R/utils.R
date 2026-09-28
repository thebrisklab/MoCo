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

make_full_data_eif <- function(a, A, Delta_Y, Delta_M, Y, gA, gDM, gDY_AXZ, eta_AXZ, eta_AXM, xi_AX, mu, pMXD, pMXZ){
  ipw_a <- as.numeric(A == a)/gA
  ipw_a_prime <- as.numeric(Delta_Y == 1&A == a)/(gDY_AXZ*gA)
  ipw_a_MD <- as.numeric(A == 0 & Delta_M == 1)/(((1 - gA)*I(a==1) + gA*I(a==0))*gDM)
  
  c_star <- pMXD/pMXZ
  
  p1 <- ipw_a_prime * c_star * (Y - mu)
  p2 <- ipw_a * (eta_AXZ - xi_AX)
  p3 <- ipw_a_MD * (eta_AXM - xi_AX)
  p4 <- xi_AX - mean(xi_AX)
  
  p1[is.na(p1)] <- 0
  p3[is.na(p3)] <- 0
  
  return(p1 + p2 + p3 + p4)
}

make_full_data_eif_easy <- function(a, A, Delta_Y, Y, gA, gDY_AXZ, xi_AX, mu) {
  ipw_a <- as.numeric(A == a) / gA
  ipw_a_prime <- as.numeric(Delta_Y == 1 & A == a) / (gDY_AXZ * gA)
  
  p1 <- ipw_a_prime * (Y - mu)
  p2 <- ipw_a * (mu - xi_AX)
  p3 <- xi_AX - mean(xi_AX)
  
  p1[is.na(p1)] <- 0
  
  return(p1 + p2 + p3)
}

#' Restore structurally missing outcome positions
#'
#' @param vec A vector or list of MoCo results after structurally missing
#'   outcome columns were removed.
#' @param seed_position Integer positions at which missing values should be
#'   restored.
#' @return An object of the same basic type as `vec`, with missing entries
#'   inserted at `seed_position`.

add_NA <- function(vec, seed_position){
  seed_position <- sort(unique(as.integer(seed_position)))
  if (length(seed_position) == 0L) return(vec)

  output_length <- length(vec) + length(seed_position)
  if (anyNA(seed_position) || any(seed_position < 1L) ||
      any(seed_position > output_length)) {
    stop("seed_position contains an invalid output position.")
  }

  keep <- setdiff(seq_len(output_length), seed_position)
  if (is.list(vec)) {
    out <- vector("list", output_length)
    out[keep] <- vec
  } else {
    out <- rep(NA, output_length)
    out[keep] <- vec
  }
  out
}

# SuperLearner screening rule that retains every supplied covariate. Some
# SuperLearner configurations request a screening function named `All`.
All <- function(Y, X, family, id = NULL, obsWeights = NULL, ...) {
  rep(TRUE, ncol(X))
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
    if (class(fit) != "error") {
      coef <- fit$solution
    } else {
      coef <- rep(0, ncol(Z))
      coef[which.min(cvRisk)] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (class(fit) != "error") {
      if (anyDupCols | anyNACols) {
        ind <- c(seq_along(coef), rmCols - 0.5)
        coef <- c(coef, rep(0, length(rmCols)))
        coef <- coef[order(ind)]
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
    if (r$status < 1 || r$status > 4) {
      warning(r$message)
    }
    if (class(r) != "error") {
      coef <- r$solution
    } else {
      coef <- rep(0, ncol(Z))
      coef[which.min(cvRisk)] <- 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] <- 0
    }
    if (anyDupCols) {
      ind <- c(seq_along(coef), dupCols - 0.5)
      coef <- c(coef, rep(0, length(dupCols)))
      coef <- coef[order(ind)]
    }
    coef[coef < 1e-04] <- 0
    coef <- coef / sum(coef)
    out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
    return(out)
  }
  list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}
