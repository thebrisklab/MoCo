# Lognormal GLM motion-density options and evaluation.

#' Fill in and check the right-hand sides of the log-normal GLM motion densities
#' @param glm_options NULL or a list with optional pMX and pMXZ right-hand-side strings.
#' @return A list containing pMX and pMXZ; each defaults to ".".
#' @noRd

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
  return(out)
}

#' Evaluate a fitted lognormal motion density
#'
#' Convert a Gaussian regression for log(M) into a density on the motion scale.
#'
#' @details The log-scale SD is sd(fit$residuals). With thresh, the supplied fit
#' must be the parent model fitted using observed motion on both sides
#' of the cutoff within the training sample. The result is f(m)/F(thresh)
#' below the cutoff and zero above it. This is not selected-sample
#' truncated-likelihood estimation. If the CDF floor activates, the
#' result is regularized and no longer integrates to one. In threshold
#' mode, empirical range_M clipping is not applied.
#'
#' @param fit Gaussian glm fit with log(M) as its response.
#' @param m Positive observed motion values at which to evaluate the density.
#' @param newdata Prediction data with one row for each element of m.
#' @param range_M NULL, or two endpoints outside which predictions are zero when thresh is NULL.
#' @param thresh NULL, or a positive cutoff for analytical conditioning on M <= thresh.
#' @param floor_cdf Lower bound for the conditioning CDF; default 0.001.
#' @return A numeric vector with one density value per element of m.
#' @noRd


lnorm_density <- function(fit, m, newdata, range_M = NULL, thresh = NULL, floor_cdf = 1e-3){
  if(length(m) == 0){
    return(numeric(0))
  }
  mean_log <- stats::predict(fit, newdata = newdata)
  sd_log <- stats::sd(fit$residuals)
  d <- (1/m) * stats::dnorm(log(m), mean = mean_log, sd = sd_log)
  if(!is.null(thresh)){
    Fc <- stats::pnorm((log(thresh) - mean_log)/sd_log)
    if(mean(Fc < floor_cdf) > 0.01){
      warning("For more than 1% of participants, the fitted P(M <= thresh) is below ", floor_cdf,
              ", so the truncated GLM density uses that floor. The motion density model may fit poorly near thresh.",
              call. = FALSE)
    }
    d <- d / pmax(Fc, floor_cdf)
    d[m > thresh] <- 0
    return(d)
  }
  if(!is.null(range_M)){
    d[m < range_M[1] | m > range_M[2]] <- 0
  }
  return(d)
}
