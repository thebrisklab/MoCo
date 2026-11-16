#' Super Learner wrapper for \code{mgcv::gam}
#'
#' A generalized additive model learner for \pkg{SuperLearner} based on
#' \code{mgcv::gam()}, used in the default \code{SL_library} of \code{moco()}.
#' Numeric covariates with more than \code{cts.num} distinct values enter
#' through penalized thin-plate regression splines \code{s()}; the other
#' covariates enter linearly. Smoothing parameters are selected by
#' \code{method} (REML by default). As for other \pkg{SuperLearner} learners,
#' if \code{mgcv::gam()} fails for a fit, \code{SuperLearner} gives
#' \code{SL.mgcv_gam} weight 0 for that fit and issues a warning.
#'
#' @param Y Outcome vector.
#' @param X Data frame of covariates for training.
#' @param newX Data frame of covariates for prediction.
#' @param family Error distribution, as passed by \code{SuperLearner}
#'   (for example \code{gaussian()} or \code{binomial()}).
#' @param obsWeights Observation weights.
#' @param cts.num Numeric covariates with more than \code{cts.num} distinct
#'   values are smoothed. Default 4.
#' @param k Maximum basis dimension of each smooth term; it is reduced for
#'   covariates with few distinct values. Default 10.
#' @param method Smoothing-parameter selection method passed to
#'   \code{mgcv::gam()}. Default \code{"REML"}.
#' @param ... Ignored.
#'
#' @return A list with elements \code{pred} (predictions for \code{newX}) and
#'   \code{fit} (an object of class \code{"SL.mgcv_gam"}), as required by
#'   \pkg{SuperLearner}.
#'
#' @examples
#' set.seed(1)
#' X <- data.frame(x1 = runif(100), x2 = rbinom(100, 1, 0.5))
#' Y <- sin(2 * pi * X$x1) + X$x2 + rnorm(100, sd = 0.3)
#' fit <- SL.mgcv_gam(Y, X, newX = X, family = gaussian(), obsWeights = rep(1, 100))
#' head(fit$pred)
#'
#' @export
SL.mgcv_gam <- function(Y, X, newX, 
                        family, obsWeights, 
                        cts.num = 4, 
                        k = 10,
                        method = "REML", ...) {
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)
  nuniq <- vapply(X, function(x) length(unique(x)), integer(1))
  cts <- vapply(X, is.numeric, logical(1)) & nuniq > cts.num
  rhs <- c(if (any(cts)) sprintf("s(%s, k = %d)", names(X)[cts], as.integer(pmin(k, nuniq[cts] - 1L))), names(X)[!cts])
  
  if (!length(rhs)) rhs <- "1"
  f <- stats::as.formula(paste("Y ~", paste(rhs, collapse = " + ")))
  fit <- mgcv::gam(f, data = data.frame(Y = Y, X), family = family,
                   weights = obsWeights, method = method)
  pred <- as.numeric(stats::predict(fit, newdata = newX, type = "response"))
  
  out <- list(object = fit)
  class(out) <- "SL.mgcv_gam"
  
  list(pred = pred, 
       fit = out)
}

#' Predictions from an \code{SL.mgcv_gam} fit
#'
#' @param object An object of class \code{"SL.mgcv_gam"}.
#' @param newdata Data frame of covariates for prediction.
#' @param ... Ignored.
#' @return A numeric vector of predictions on the response scale.
#' @export
predict.SL.mgcv_gam <- function(object, newdata, ...) {
  as.numeric(stats::predict(object$object, newdata = as.data.frame(newdata), type = "response"))
}
