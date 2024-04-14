#' Evaluate EIF for stochastic interventional (in)direct effects
#'
#' @param a Either 0 or 1.
#' @param A A binary vector of length n \code{number of participants} representing the diagnosis status.
#' @param Delta A binary vector indicating data usability, such as whether data pass a motion criteria.
#' @param Y A vector of continuous outcome of interest.
#' @param gA Estimate of P(A = 1 | X_i), i = 1,...,n.
#' @param gD Estimate of P(Delta = 1 | A_i, X_i), i = 1,...,n.
#' @param eta_AXZ Estimate of E(mu_pseudo_A | A_i, Delta_i, X_i).
#' @param eta_AXM Estimate of E(mu_pseudo_A_star | A_i, M_i, X_i).
#' @param xi_AX Estimate of E(eta_AXZ | A_i, Delta_i, X_i).
#' @param mu Estimate of E(Y| A_i, M_i, X_i, Z_i).
#' @param pMX Estimate of p(M | 0, Delta_i = 1, X_i).
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

#' Add NA value back indicating seed position
#'
#' @param vec A vector indicating the MoCo estimates, adjusted association, z-scores, or significant regions.
#' @param seed_position Position indicator of where the seed region is.
#' @return A vector with NA values inserted to indicate the position of the seed region.

add_NA <- function(vec, seed_position){
  if(seed_position == (length(vec) + 1)){
    vec <- c(vec[1:(seed_position-1)], NA)
  }else{
    vec <- c(vec[1:(seed_position-1)], NA, vec[seed_position:length(vec)])
  }
  
  return(vec)
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

