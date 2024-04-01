#' Evaluate EIF for stochastic interventional (in)direct effects
#' 
#' @param a either 0 or 1
#' @param A a binary vector of length n \code{number of participants}, denoting diagnosis status
#' @param Delta a binary vector of data usability. For example, this corresponds to the binary variable indicating 
#' whether data pass a motion criteria, such that these data would be excluded in conventional analyses.
#' @param Y a vector of continuous outcome of interest
#' 
#' @param gA estimate of P(A = 1 | X_i), i = 1,...,n
#' @param gD estimate of P(Delta = 1 | A_i, X_i), i = 1,...,n
#' @param eta_AXZ estimate of E(mu_pseudo_A | A_i, Delta_i, X_i)
#' @param eta_AXM estimate of E(mu_pseudo_A_star | A_i, M_i, X_i)
#' @param xi_AX estimate of E(eta_AXZ | A_i, Delta_i, X_i)
#' @param mu estimate of E(Y| A_i, M_i, X_i, Z_i)
#' @param pMX estimate of p(M | 0, Delta_i = 1, X_i)
#' @param pMXZ Estimate of p(M | A_i, X_i, Z_i)
#' 
#' @return An n-length vector of the estimated EIF evaluated on the observations

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


#' Define a simulation setting where some nuisance quantities are consistently estimated
#' 
#' If a nuisance quantity is included in \code{which_correct}, the proper specification is 
#' given for the regression function. If a nuisance quantity is not included in \code{which_correct},
#' an intercept-only model is used instead. 
#' 
#' @param which_correct A character vector indicating which nuisance regressions are consistently estimated.
#' @return A named list that can be passed to \code{one_sim} function to run a simulation 
#' under a certain pattern of nuisance parameter (mis)specification.

create_setting = function(which_correct){
  # propensity score: disease status, conditional on baseline covariates
  glm_gA = ifelse("gA" %in% which_correct, "X", 1)
  # binary mediator regression
  glm_gDM = ifelse("gDM" %in% which_correct, "A*X", 1)

  # missingness indicator regression
  # if("gDY_X" %in% which_correct){glm_gDY_X = NULL; SL_gDY_X = SL_library}else{glm_gDY_X = 1; SL_gDY_X = NULL}
  # if("gDY_AX" %in% which_correct){glm_gDY_AX = NULL; SL_gDY_AX = SL_library}else{glm_gDY_AX = 1; SL_gDY_AX = NULL}
  # if("gDY_AXZ" %in% which_correct){glm_gDY_AXZ = NULL; SL_gDY_AXZ = SL_library}else{glm_gDY_AXZ = 1; SL_gDY_AXZ = NULL}
  if("gDY_AX" %in% which_correct){glm_gDY_AX = "A*X"}else{glm_gDY_AX = 1}
  if("gDY_AXZ" %in% which_correct){glm_gDY_AXZ = "A*X*Z"}else{glm_gDY_AXZ = 1}
  
  # outcome regression 
  glm_mu_AMXZ = ifelse("mu_AMXZ" %in% which_correct, "A + M + X + Z", 1)
  
  # density estimation for mediator, conditioning on disease status, binary mediator, and baseline covariates
  if("pMX" %in% which_correct){HAL_pMX = TRUE; glm_pMX = NULL}else{HAL_pMX = FALSE; glm_pMX = 1}
  # density estimation for mediator, conditioning on disease status, binary mediator, baseline covariates, and mediator-outcome confounder
  if("pMXZ" %in% which_correct){HAL_pMXZ = TRUE; glm_pMXZ = NULL}else{HAL_pMXZ = FALSE; glm_pMXZ = 1}
  
  # fit pseudo-outcome regression, conditioning on disease status, mediator and baseline covariates
  # if("eta_AXZ" %in% which_correct){glm_eta_AXZ = NULL;SL_eta_AXZ = SL_library}else{glm_eta_AXZ = 1;SL_eta_AXZ = NULL}
  # if("xi_AX" %in% which_correct){glm_xi_AX = NULL;SL_xi_AX = SL_library}else{glm_xi_AX = 1;SL_xi_AX = NULL}
  # if("eta_AXM" %in% which_correct){glm_eta_AXM = NULL; SL_eta_AXM = NULL}else{glm_eta_AXM = 1; SL_eta_AXM = NULL} # HAL
  if("eta_AXZ" %in% which_correct){glm_eta_AXZ = "A*X*Z"}else{glm_eta_AXZ = 1}
  if("xi_AX" %in% which_correct){glm_xi_AX = "A*X"}else{glm_xi_AX = 1}
  if("eta_AXM" %in% which_correct){glm_eta_AXM = "A*M*X"}else{glm_eta_AXM = 1} 
  
  return(
    list(glm_formula = list(gA = glm_gA,
                            gDM = glm_gDM,
                            gDY_AX = glm_gDY_AX,
                            gDY_AXZ = glm_gDY_AXZ,
                            mu_AMXZ = glm_mu_AMXZ,
                            eta_AXZ = glm_eta_AXZ,
                            xi_AX = glm_xi_AX,
                            eta_AXM = glm_eta_AXM,
                            pMX = glm_pMX,
                            pMXZ = glm_pMXZ),
         HAL_pMX = HAL_pMX,
         HAL_pMXZ = HAL_pMXZ
    ))
}

#' Temporary fix for convex combination method mean squared error
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another
tmp_method.CC_LS = function() {
  computeCoef = function(Z, Y, libraryNames, verbose, obsWeights,
                          errorsInLibrary = NULL, ...) {
    cvRisk = apply(Z, 2, function(x) {
      mean(obsWeights * (x -
                           Y)^2)
    })
    names(cvRisk) = libraryNames
    compute = function(x, y, wt = rep(1, length(y))) {
      wX = sqrt(wt) * x
      wY = sqrt(wt) * y
      D = crossprod(wX)
      d = crossprod(wX, wY)
      A = cbind(rep(1, ncol(wX)), diag(ncol(wX)))
      bvec = c(1, rep(0, ncol(wX)))
      fit = tryCatch(
        {
          quadprog::solve.QP(
            Dmat = D, dvec = d, Amat = A,
            bvec = bvec, meq = 1
          )
        },
        error = function(e) {
          out = list()
          class(out) = "error"
          out
        }
      )
      invisible(fit)
    }
    modZ = Z
    naCols = which(apply(Z, 2, function(z) {
      all(z == 0)
    }))
    anyNACols = length(naCols) > 0
    if (anyNACols) {
      warning(paste0(
        paste0(libraryNames[naCols], collapse = ", "),
        " have NAs.", "Removing from super learner."
      ))
    }
    tol = 4
    dupCols = which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols = length(dupCols) > 0
    # if (anyDupCols) {
    #   warning(paste0(
    #     paste0(libraryNames[dupCols], collapse = ", "),
    #     " are duplicates of previous learners.", " Removing from super learner."
    #   ))
    # }
    if (anyDupCols | anyNACols) {
      rmCols = unique(c(naCols, dupCols))
      modZ = Z[, -rmCols, drop = FALSE]
    }
    fit = compute(x = modZ, y = Y, wt = obsWeights)
    if (class(fit) != "error") {
      coef = fit$solution
    } else {
      coef = rep(0, ncol(Z))
      coef[which.min(cvRisk)] = 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] = 0
    }
    if (class(fit) != "error") {
      if (anyDupCols | anyNACols) {
        ind = c(seq_along(coef), rmCols - 0.5)
        coef = c(coef, rep(0, length(rmCols)))
        coef = coef[order(ind)]
      }
      coef[coef < 1e-04] = 0
      coef = coef / sum(coef)
    }
    if (!sum(coef) > 0) {
      warning("All algorithms have zero weight", call. = FALSE)
    }
    list(cvRisk = cvRisk, coef = coef, optimizer = fit)
  }
  computePred = function(predY, coef, ...) {
    predY %*% matrix(coef)
  }
  out = list(
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
tmp_method.CC_nloglik = function() {
  computePred = function(predY, coef, control, ...) {
    if (sum(coef != 0) == 0) {
      stop("All metalearner coefficients are zero, cannot compute prediction.")
    }
    stats::plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*%
                    matrix(coef[coef != 0]))
  }
  computeCoef = function(Z, Y, libraryNames, obsWeights, control,
                          verbose, ...) {
    tol = 4
    dupCols = which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols = length(dupCols) > 0
    modZ = Z
    if (anyDupCols) {
      # warning(paste0(
      #   paste0(libraryNames[dupCols], collapse = ", "),
      #   " are duplicates of previous learners.", " Removing from super learner."
      # ))
      modZ = modZ[, -dupCols, drop = FALSE]
    }
    modlogitZ = trimLogit(modZ, control$trimLogit)
    logitZ = trimLogit(Z, control$trimLogit)
    cvRisk = apply(logitZ, 2, function(x) {
      -sum(2 * obsWeights *
             ifelse(Y, stats::plogis(x, log.p = TRUE), stats::plogis(x,
                                                                     log.p = TRUE,
                                                                     lower.tail = FALSE
             )))
    })
    names(cvRisk) = libraryNames
    obj_and_grad = function(y, x, w = NULL) {
      y = y
      x = x
      function(beta) {
        xB = x %*% cbind(beta)
        loglik = y * stats::plogis(xB, log.p = TRUE) + (1 -
                                                           y) * stats::plogis(xB, log.p = TRUE, lower.tail = FALSE)
        if (!is.null(w)) {
          loglik = loglik * w
        }
        obj = -2 * sum(loglik)
        p = stats::plogis(xB)
        grad = if (is.null(w)) {
          2 * crossprod(x, cbind(p - y))
        } else {
          2 * crossprod(x, w * cbind(p - y))
        }
        list(objective = obj, gradient = grad)
      }
    }
    lower_bounds = rep(0, ncol(modZ))
    upper_bounds = rep(1, ncol(modZ))
    if (anyNA(cvRisk)) {
      upper_bounds[is.na(cvRisk)] = 0
    }
    r = tryCatch(
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
        out = list()
        class(out) = "error"
        out
      }
    )
    if (r$status < 1 || r$status > 4) {
      warning(r$message)
    }
    if (class(r) != "error") {
      coef = r$solution
    } else {
      coef = rep(0, ncol(Z))
      coef[which.min(cvRisk)] = 1
    }
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] = 0
    }
    if (anyDupCols) {
      ind = c(seq_along(coef), dupCols - 0.5)
      coef = c(coef, rep(0, length(dupCols)))
      coef = coef[order(ind)]
    }
    coef[coef < 1e-04] = 0
    coef = coef / sum(coef)
    out = list(cvRisk = cvRisk, coef = coef, optimizer = r)
    return(out)
  }
  list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}

