# Small shared helpers for MoCo.

# ---------------------------------------------------------------------------
# output formatting
# ---------------------------------------------------------------------------

#' Add NA value back indicating seed position
#'
#' @param vec A vector indicating the MoCo estimates, adjusted association, z-scores, or significant regions.
#' @param seed_position Position indicator of where the seed region is.
#' @return A vector with NA values inserted to indicate the position of the seed region.

add_NA <- function(vec, seed_position){
  # insert NA at each seed position (positions refer to the original column indices of Y)
  for(pos in sort(seed_position)){
    vec <- append(vec, NA, after = pos - 1)
  }
  
  return(vec)
}

# ---------------------------------------------------------------------------
# motion density helpers (shared by the glm, gamlss and hal density methods)
# ---------------------------------------------------------------------------

# resolve density_method into one method for pMX/pMXD and one for pMXZ/pMXZD
resolve_density_method <- function(density_method){
  if(is.null(density_method)){
    density_method <- "hal"
  }
  nms <- names(density_method)
  if(!is.character(density_method) || !(length(density_method) %in% c(1, 2)) ||
     !all(tolower(density_method) %in% c("glm", "hal", "gamlss"))){
    stop("density_method must be \"hal\", \"glm\" or \"gamlss\", or a length-2 vector c(pMX = , pMXZ = ) of these.")
  }
  density_method <- tolower(density_method)
  names(density_method) <- nms
  if(length(density_method) == 1){
    density_method <- rep(density_method, 2)
  }
  if(!is.null(names(density_method)) && all(names(density_method) %in% c("pMX", "pMXZ"))){
    density_method <- density_method[c("pMX", "pMXZ")]
  }
  names(density_method) <- c("pMX", "pMXZ")
  density_method
}

# fill in and check the right-hand sides of the log-normal GLM motion densities
check_glm_options <- function(glm_options){
  out <- list(pMX = ".", pMXZ = ".")
  if(is.null(glm_options)){
    return(out)
  }
  if(!is.list(glm_options) || (length(glm_options) > 0 && is.null(names(glm_options))) ||
     !all(names(glm_options) %in% c("pMX", "pMXZ"))){
    stop("glm_options must be a list with entries pMX and/or pMXZ, e.g. list(pMX = \".\", pMXZ = \".\").")
  }
  for(nm in names(glm_options)){
    if(!is.null(glm_options[[nm]])){
      if(!is.character(glm_options[[nm]]) || length(glm_options[[nm]]) != 1){
        stop("glm_options$", nm, " must be one character string, the right-hand side of the formula (e.g. \".\").")
      }
      out[[nm]] <- glm_options[[nm]]
    }
  }
  out
}

# summary of the density-ratio weights pMXD / pMXZ that enter the efficient influence function
weight_summary <- function(w){
  w <- w[is.finite(w)]
  if(length(w) == 0){
    return(c(n = 0, max_weight = NA, ess_ratio = NA))
  }
  c(n = length(w), max_weight = max(w), ess_ratio = sum(w)^2 / (length(w) * sum(w^2)))
}

# ---------------------------------------------------------------------------
# SuperLearner metalearner fixes
# ---------------------------------------------------------------------------

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
    if (!inherits(r, "error") && (r$status < 1 || r$status > 4)) {
      warning(r$message)
    }
    if (!inherits(r, "error")) {
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
