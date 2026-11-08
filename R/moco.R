#' Non-parametric efficient motion-controlled brain phenotype estimators
#' 
#' @details Compute non-parametric efficient motion-controlled brain phenotype estimators (MoCo)
#' 
#' @param X A dataframe or matrix containing demographic confounders that would ideally be balanced in a randomized controlled trial.
#' @param Z A dataframe or matrix of covariates representing brain phenotypes.
#' @param A A binary vector of length n (number of participants), serving as a group indicator, such as diagnosis group or control group.
#' @param M A numeric vector of length n representing continuous motion values for each participant.
#' @param Y A matrix of dimension n x p, where n is the number of participants, and p is the number of regions of interest.
#'          If it represents seed-based association measures:
#'            Each (i, j) element denotes participant i's association measure between the seed region and region j.
#'            The column representing the association measure of the seed region with itself should be filled with NA values to indicate its position.
#'          If it represents other types of association measures:
#'            Each (i, j) element denotes participant i's association measure between two brain regions of interest, such as the upper diagonal part of the functional connectivity matrix.
#'            No NA values are allowed in \code{Y} in this case.
#' @param Delta_M Required. A binary (0/1) vector of length n indicating whether motion is available and meets the motion quality control (QC) inclusion criteria.
#'          Set \code{Delta_M} = 1 if the participant passes motion QC; otherwise, set \code{Delta_M} = 0.
#'          Any QC definition can be used, e.g., a threshold on mean framewise displacement (\code{Delta_M = as.numeric(M < 0.2)}),
#'          or criteria based on scrubbing or the number of remaining frames. \code{Delta_M} must be 0 when \code{M} is missing.
#' @param Delta_Y A binary vector indicating the non-missingness and whether the brain image data \code{Y} passes quality control after preprocessing. 
#'          Set \code{Delta_Y = 1} if \code{Y} is usable; otherwise, set \code{Delta_Y = 0}.
#' 
#' @param SL_library SuperLearner library for estimating nuisance regressions. 
#'                   Defaults to c("SL.earth","SL.glmnet","SL.gam","SL.glm", "SL.glm.interaction", "SL.step","SL.step.interaction","SL.xgboost","SL.ranger","SL.mean") if not specified.
#' @param SL_library_customize Customize SuperLearner library for estimating each nuisance regression. 
#'                    - \code{gA}: SuperLearner library for estimating the propensity score.
#'                    - \code{gDM}: SuperLearner library for estimating the probability P(Delta_M = 1 | A, X).
#'                    - \code{gDY_AX}: SuperLearner library for estimating the probability P(Delta_Y = 1 | A, X).
#'                    - \code{gDY_AXZ}: SuperLearner library for estimating the probability P(Delta_Y = 1 | A, X, Z).
#'                    - \code{mu_AMXZ}: SuperLearner library for estimating the outcome regression E(Y | Delta_Y = 1, A, M, X, Z).
#'                    - \code{eta_AXZ}: SuperLearner library for estimating E(mu_AMXZ pMXD / pMXZD | A, X, Z, Delta_M = 1).                
#'                    - \code{eta_AXM}: SuperLearner library for estimating E(mu_AMXZ pMX/pMXZ gDY_AX/gDY_AXZ | A, M, X, Delta_Y = 1).
#'                    - \code{xi_AX}: SuperLearner library for estimating E(eta_AXZ | A, X).
#'                    
#' @param glm_formula All glm formulas default to NULL, indicating SuperLearner will be used for nuisance regressions.
#'                    - \code{gA}: GLM formula for estimating the propensity score.
#'                    - \code{gDM}: GLM formula for estimating the probability P(Delta_M = 1 | A, X).
#'                    - \code{gDY_AX}: GLM formula for estimating the probability P(Delta_Y = 1 | A, X).
#'                    - \code{gDY_AXZ}: GLM formula for estimating the probability P(Delta_Y = 1 | A, X, Z).
#'                    - \code{mu_AMXZ}: GLM formula for estimating the outcome regression E(Y | Delta_Y = 1, A, M, X, Z).
#'                    - \code{eta_AXZ}: GLM formula for estimating E(mu_AMXZ pMXD / pMXZD | A, X, Z, Delta_M = 1).                
#'                    - \code{eta_AXM}: GLM formula for estimating E(mu_AMXZ pMX/pMXZ gDY_AX/gDY_AXZ | A, M, X, Delta_Y = 1).
#'                    - \code{xi_AX}: GLM formula for estimating E(eta_AXZ | A, X).
#'                    
#' @param density_method Method for the motion densities p(m | a, x, Delta) and p(m | a, x, z, Delta), with Delta = Delta_Y or Delta_M:
#'          \code{"hal"} (default; highly adaptive lasso, \code{haldensify}, controlled by \code{HAL_options}),
#'          \code{"glm"} (log-normal regression, controlled by \code{glm_options}), or
#'          \code{"gamlss"} (generalized gamma GAMLSS, controlled by \code{gamlss_options} and \code{Delta_M_threshold}).
#'          Upper or lower case is accepted. A length-2 vector \code{c(pMX = , pMXZ = )} sets the two density pairs separately.
#'          Options for a method that is not used are ignored, with a message. The former arguments \code{HAL_pMX} and
#'          \code{HAL_pMXZ} are deprecated: when supplied, they are translated to \code{density_method} with a warning.
#' @param HAL_options Options for \code{density_method = "hal"}, passed to \code{haldensify::haldensify()}:
#'                   - \code{max_degree}: The highest order of interaction terms for generating basis functions.
#'                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression.
#'                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions.
#' @param glm_options Right-hand sides of the log-normal GLM motion densities, used only when \code{density_method = "glm"}:
#'                   - \code{pMX}: formula for p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M = 1). Default \code{"."} (all of \code{A} and \code{X}).
#'                   - \code{pMXZ}: formula for p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M = 1). Default \code{"."} (all of \code{A}, \code{X} and \code{Z}).
#'          The former entries \code{glm_formula$pMX} and \code{glm_formula$pMXZ} are deprecated: when supplied, they are moved to \code{glm_options} with a warning.
#' @param Delta_M_threshold Used only with \code{density_method = "gamlss"}. When \code{Delta_M = as.numeric(M < c)}, set \code{Delta_M_threshold = c}:
#'          p(m | a, x, Delta_M = 1) and p(m | a, x, z, Delta_M = 1) are then fitted with a generalized gamma right-truncated at c (see \code{\link{GGtr}}),
#'          and set to 0 for M >= c. Defaults to \code{NULL} (Delta_M is not a cut on M; an ordinary generalized gamma is fitted on the Delta_M = 1 participants
#'          and set to 0 outside their range of M). With \code{NULL}, MoCo warns when Delta_M looks like a threshold on M.
#' @param gamlss_options Options for \code{density_method = "gamlss"} (fitted with \code{gamlss} and \code{gamlss.dist}, which are installed with MoCo). Entries left out take their defaults:
#'                   - \code{family}: \code{"GG"} (generalized gamma; the only family currently supported).
#'                   - \code{mu}: \code{"linear"} (all covariates entered linearly, default), \code{"smooth"} (continuous covariates entered as \code{pb()} P-splines, the others linearly), or a one-sided formula.
#'                   - \code{sigma}: \code{"linear"} (default), \code{"smooth"}, \code{"1"} (constant), \code{"A"} (group only), or a one-sided formula.
#'                   - \code{nu}: \code{"1"} (constant nu, estimated; default) or \code{0} (nu fixed at 0, i.e. the log-normal distribution).
#'                   - \code{smooth}: names of the continuous covariates that get \code{pb()} when \code{mu} or \code{sigma} is \code{"smooth"}.
#'                     \code{NULL} (default) uses every numeric covariate with more than 10 distinct values.
#'                   - \code{formula}: optional \code{list(pMX =, pMXD =, pMXZ =, pMXZD =)}, each a \code{list(mu =, sigma =, nu =)} of one-sided formulas (nu may be \code{0}); overrides the settings above for that density.
#'                   - \code{method}: algorithms tried in order until one converges, any of \code{"RS"}, \code{"CG"}, \code{"mixed"}. Default \code{c("RS", "CG", "mixed")}.
#'                   - \code{n.cyc}, \code{c.crit}: maximum number of outer cycles (default 500) and convergence criterion (default 0.001), passed to \code{gamlss.control}.
#'                   - \code{nu0_fallback}: if no algorithm converges with nu estimated, refit with nu fixed at 0. Default \code{TRUE}. If that also fails, the GLM log-normal density is used for that component, with a warning.
#'                   - \code{clamp}: clamp continuous covariates to their training range before prediction. Default \code{FALSE}.
#'                   - \code{trace}: print the gamlss iterations. Default \code{FALSE}.
#' 
#' @param cv_folds Number of folds used for cross-fitting when \code{cross_fit = TRUE}. Defaults to 5.
#' @param cross_fit Logical indicating whether to develop the estimator based on cross-fitting. Defaults to TRUE.
#' @param test Logical indicating whether to conduct hypothesis testing based on simultaneous confidence band. Defaults to TRUE.
#' @param fwer A vector of family-wise error rates (FWER) to control for multiple hypothesis testing. Defaults to c(0.05). Set to NULL if test is FALSE.
#' @param hypo_method Critical-value simulation: "mvn" (default) or "multiplier".
#' @param hypo_n_sim Number of Monte Carlo draws for the simultaneous test. Defaults to 100000, as in the original MoCo; it can be lowered (e.g. to 10000) for quick exploratory runs with many outcomes.
#' @param hypo_chunk_size Draws processed per simulation batch.
#' @param hypo_seed Critical-value simulation seed, separate from nuisance seeds.
#' @param seed_rng Specifies the value of seed(s) for nuisance regression calculation using super learner. Can be a vector. Defaults to value 1. When several seeds are given, the full analysis is run once per seed; the estimates and the participant-level efficient influence functions are averaged across seeds, the covariance is recomputed from the averaged efficient influence functions, and one simultaneous test is conducted across all outcomes in \code{Y}.
#' @param ... Used only for the deprecated \code{HAL_pMX} and \code{HAL_pMXZ} (see \code{density_method}). Supplying the removed arguments \code{thresh} or \code{seed_rgn} gives an informative error.
#' 
#' @import SuperLearner 
#' @import haldensify 
#' @import MASS
#' @importFrom stats binomial gaussian dnorm sd cov
#' 
#' @return A list with named entries 
#' \describe{
#'   \item{est}{A two times p matrix showing the one-step estimators of the control group and the disease group for each functional connectivity of interest, respectively.}
#'   \item{adj_association}{A p-length vector showing the motion-adjusted association for each functional connectivity of interest, respectively.}
#'   \item{z_score}{A p-length vector of z_score for each functional connectivity if \code{test} is TRUE.}
#'   \item{fwer_critical_value}{Critical value(s) of the simultaneous confidence band (the (1 - FWER) quantile of the maximum absolute z-statistic) for each specified FWER.}
#'   \item{significant_regions}{A p-length vector of TRUE or FALSE, indicating significance for each functional connectivity if \code{test} is TRUE.}
#'   \item{density_log}{With \code{density_method = "gamlss"}: one row per seed, fold and density with the algorithm used, whether nu was fixed at 0, truncation, convergence, iterations, the fitted nu, the sample size, and whether a fallback was needed. \code{NULL} otherwise.}
#'   \item{weight_diagnostics}{One row per seed and group (A0, A1) with the number, maximum and effective-sample-size ratio of the density-ratio weights pMXD / pMXZ in the efficient influence function. MoCo warns when an effective-sample-size ratio is below 0.1.}
#' }
#' 
#' @examples
#' data(moco_example)
#' d <- moco_example
#'
#' # Fast example (a few seconds): one outcome (Frontoparietal) plus the
#' # seed-to-self column, a small Super Learner library and GAMLSS motion densities
#' fit <- moco(X = d$X, Z = d$Z, A = d$A, M = d$M, Y = d$Y[, 6:7],
#'             Delta_M = d$Delta_M, Delta_Y = d$Delta_Y,
#'             SL_library = c("SL.mean", "SL.glm"),
#'             density_method = "gamlss", cv_folds = 5, seed_rng = 1)
#' fit$adj_association
#' fit$significant_regions
#' fit$density_log
#'
#' # The full GLM, GAMLSS and HAL analysis of the README tutorial (about 50 min
#' # in total) is in an installed script:
#' # file.show(system.file("examples", "full_analysis.R", package = "MoCo"))
#'
#' @export

moco <- function(
  X, Z, A, M, Y, 
  Delta_M, 
  Delta_Y,
  SL_library = c("SL.earth","SL.glmnet","SL.gam","SL.glm", "SL.glm.interaction", "SL.step","SL.step.interaction","SL.xgboost","SL.ranger","SL.mean"),
  SL_library_customize = list(
    gA = NULL, 
    gDM = NULL,
    gDY_AX = NULL,
    gDY_AXZ = NULL,
    mu_AMXZ = NULL,
    eta_AXZ = NULL,
    eta_AXM = NULL,
    xi_AX = NULL
  ), 
  glm_formula = list(gA = NULL, 
                     gDM = NULL,
                     gDY_AX = NULL,
                     gDY_AXZ = NULL,
                     mu_AMXZ = NULL,
                     eta_AXZ = NULL,
                     eta_AXM = NULL,
                     xi_AX = NULL),
  density_method = "hal",
  HAL_options = list(max_degree = 3, lambda_seq = exp(seq(-1, -10, length = 100)), num_knots = c(1000, 500, 250)),
  glm_options = list(pMX = ".", pMXZ = "."),
  Delta_M_threshold = NULL,
  gamlss_options = NULL,
  cross_fit = TRUE,
  cv_folds = 5,
  test = TRUE,
  fwer = 0.05, 
  seed_rng = 1, 
  hypo_method = "mvn",
  hypo_n_sim = 100000L,
  hypo_chunk_size = 1000L,
  hypo_seed = 1,
  ...
){
  
  # which density options the user supplied (recorded before any of them is changed)
  density_method_given <- !missing(density_method)
  HAL_options_given <- !missing(HAL_options)
  glm_options_given <- !missing(glm_options)
  
  # arguments that were renamed or removed in this version
  dot_args <- list(...)
  dots <- names(dot_args)
  if("thresh" %in% dots){
    stop("The argument 'thresh' has been removed. Please define Delta_M directly, e.g., Delta_M = as.numeric(M < 0.2).")
  }
  if("seed_rgn" %in% dots){
    stop("The argument 'seed_rgn' has been renamed to 'seed_rng'.")
  }
  
  # deprecated HAL_pMX / HAL_pMXZ: translated to density_method
  if(any(c("HAL_pMX", "HAL_pMXZ") %in% dots)){
    if(!density_method_given){
      hal_pMX <- if(is.null(dot_args[["HAL_pMX"]])) TRUE else isTRUE(dot_args[["HAL_pMX"]])
      hal_pMXZ <- if(is.null(dot_args[["HAL_pMXZ"]])) TRUE else isTRUE(dot_args[["HAL_pMXZ"]])
      density_method <- c(pMX = if(hal_pMX) "hal" else "glm", pMXZ = if(hal_pMXZ) "hal" else "glm")
      warning("HAL_pMX and HAL_pMXZ are deprecated; please use density_method. They were translated to density_method = c(pMX = \"",
              density_method[["pMX"]], "\", pMXZ = \"", density_method[["pMXZ"]], "\").", call. = FALSE)
      density_method_given <- TRUE
    }else{
      warning("HAL_pMX and HAL_pMXZ are deprecated and are ignored because density_method is given.", call. = FALSE)
    }
  }
  
  # deprecated glm_formula$pMX / glm_formula$pMXZ: moved to glm_options
  if(!is.null(glm_formula[["pMX"]]) || !is.null(glm_formula[["pMXZ"]])){
    if(!glm_options_given){
      glm_options <- list(pMX = glm_formula[["pMX"]], pMXZ = glm_formula[["pMXZ"]])
      warning("glm_formula$pMX and glm_formula$pMXZ are deprecated; please use glm_options. They were moved to glm_options.", call. = FALSE)
      glm_options_given <- TRUE
    }else{
      warning("glm_formula$pMX and glm_formula$pMXZ are deprecated and are ignored because glm_options is given.", call. = FALSE)
    }
  }
  glm_formula[["pMX"]] <- NULL
  glm_formula[["pMXZ"]] <- NULL
  glm_options <- check_glm_options(glm_options)
  
  # total number of runs for handling Monte Carlo variability
  r <- length(seed_rng)
  
  # check if A, X, Z are complete
  if(sum(is.na(A)) != 0 | sum(is.na(X)) != 0 | sum(is.na(Z)) != 0) {
    stop("There are missing values in A, X, or Z. Please input full data for these variables.")
  }
  
  # change Y to matrix if it is a vector
  if(is.null(dim(Y))){
    Y <- as.matrix(Y)
  }
  n <- nrow(Y)
  
  # check Delta_M: required for every definition of motion quality control
  if(missing(Delta_M) || is.null(Delta_M)){
    stop("Delta_M is required. Please provide a 0/1 vector indicating whether each participant passes motion quality control.")
  }
  if(length(Delta_M) != n){
    stop("Delta_M must have length n = nrow(Y).")
  }
  if(any(is.na(Delta_M)) || !all(Delta_M %in% c(0, 1))){
    stop("Delta_M must be a 0/1 vector without missing values.")
  }
  
  # check Delta_Y
  if(length(Delta_Y) != n){
    stop("Delta_Y must have length n = nrow(Y).")
  }
  if(any(is.na(Delta_Y)) || !all(Delta_Y %in% c(0, 1))){
    stop("Delta_Y must be a 0/1 vector without missing values.")
  }
  
  # check M: missing motion is only allowed for participants with Delta_M = 0 and Delta_Y = 0
  if(!is.null(M)){
    if(length(M) != n){
      stop("M must have length n = nrow(Y).")
    }
    if(any(is.na(M) & Delta_M == 1)){
      stop("Delta_M = 1 for participants with missing M. Delta_M must be 0 when M is missing.")
    }
    if(any(is.na(M) & Delta_Y == 1)){
      stop("Delta_Y = 1 for participants with missing M. Participants with usable Y must have observed M.")
    }
  }
  
  # motion density method and the threshold that defines Delta_M
  dens_method <- resolve_density_method(density_method)
  if(!is.null(M)){
    if(HAL_options_given && !any(dens_method == "hal")){
      message("HAL_options is only used by density_method = \"hal\"; it is ignored here.")
    }
    if(glm_options_given && !any(dens_method == "glm")){
      message("glm_options is only used by density_method = \"glm\"; it is ignored here.")
    }
    if(!is.null(gamlss_options) && !any(dens_method == "gamlss")){
      message("gamlss_options is only used by density_method = \"gamlss\"; it is ignored here.")
    }
  }
  uses_gamlss <- !is.null(M) && any(dens_method == "gamlss")
  if(uses_gamlss){
    gg_check_packages()
    gg_check_options(gamlss_options)
  }
  Delta_M_threshold <- check_Delta_M_threshold(Delta_M_threshold, M, Delta_M, uses_gamlss)
  
  # seed position
  seed_position <- which(apply(Y, 2, function(x){
    sum(is.na(x))
  }) == nrow(Y))
  
  if(length(seed_position) != 0){
    # remove NA column for subsequent calculation
    Y <- Y[, -seed_position, drop = FALSE]
  }
  
  # check Y: participants with Delta_Y = 1 must have observed Y
  if(any(is.na(Y[Delta_Y == 1, , drop = FALSE]))){
    stop("Y has missing values for participants with Delta_Y = 1. Set Delta_Y = 0 for participants with unusable Y.")
  }
  
  # create matrix and vector for storing fitting results
  est <- list()
  adj_association <- matrix(nrow = r, ncol = ncol(Y))
  
  # running sum of the participant-level efficient influence functions across runs
  eif_sum <- NULL
  
  # density fit log and density-ratio weight diagnostics across runs
  density_log <- NULL
  weight_diagnostics <- NULL
  
  for(i in 1:r){
    
    # whether to use cross-fit or not
    if(cross_fit){
      result <- one_step_cross(
        X = X, Z = Z, A = A, M = M, Y = Y, 
        Delta_M = Delta_M, 
        Delta_Y = Delta_Y,
        SL_library = SL_library,
        SL_library_customize = SL_library_customize,
        glm_formula = glm_formula,
        HAL_options = HAL_options,
        glm_options = glm_options,
        density_method = dens_method,
        Delta_M_threshold = Delta_M_threshold,
        gamlss_options = gamlss_options,
        seed = seed_rng[i], 
        cv_folds = cv_folds
      )
    }else{
      result <- one_step(
        X = X, Z = Z, A = A, M = M, Y = Y, 
        Delta_M = Delta_M, 
        Delta_Y = Delta_Y,
        SL_library = SL_library,
        SL_library_customize = SL_library_customize,
        glm_formula = glm_formula,
        HAL_options = HAL_options,
        glm_options = glm_options,
        density_method = dens_method,
        Delta_M_threshold = Delta_M_threshold,
        gamlss_options = gamlss_options,
        seed = seed_rng[i]
      )
    }
    
    # store estimation results
    est <- c(est, list(result$est))
    adj_association[i,] <- result$adj_association
    if(!is.null(result$density_log)){
      density_log <- rbind(density_log, data.frame(seed = seed_rng[i], result$density_log, row.names = NULL))
    }
    if(!is.null(result$weight_diagnostics)){
      weight_diagnostics <- rbind(weight_diagnostics, data.frame(seed = seed_rng[i], group = rownames(result$weight_diagnostics),
                                                                 result$weight_diagnostics, row.names = NULL))
    }
    
    # accumulate participant-level efficient influence functions for hypothesis testing
    if(test){
      if(is.null(eif_sum)){
        eif_sum <- result$eif_mat
      }else{
        eif_sum <- Map(`+`, eif_sum, result$eif_mat)
      }
    }
  }
  
  # flag unstable density-ratio weights
  if(!is.null(weight_diagnostics) && any(weight_diagnostics$ess_ratio < 0.1, na.rm = TRUE)){
    warning("The effective-sample-size ratio of the density-ratio weights pMXD / pMXZ is below 0.1 for at least one group and seed ",
            "(see weight_diagnostics). The motion density models may be too flexible or misspecified.", call. = FALSE)
  }
  
  # merge multi-run estimation results
  est <- Reduce("+", est)/r
  if(length(seed_position) != 0){
    est <- t(apply(est, 1, add_NA, seed_position))
  }
  est <- as.data.frame(est, row.names = c('est_A0', 'est_A1'))
  colnames(est) <- NULL
  
  # motion-controlled association
  adj_association <- colMeans(adj_association)
  
  if(test){
    # average the participant-level efficient influence functions across runs,
    # recompute their covariance, and conduct one simultaneous test across all outcomes
    eif_avg <- lapply(eif_sum, function(x){x / r})
    cov_avg <- lapply(eif_avg, function(x){stats::cov(x) / nrow(x)})
    hypo <- hypo_test(
      list(adj_association = adj_association, eif_mat = eif_avg, cov_mat = cov_avg), 
      fwer = fwer, seed = hypo_seed, method = hypo_method,
      n_sim = hypo_n_sim, chunk_size = hypo_chunk_size
    )
    
    z_score <- hypo$z_score
    if(length(seed_position) != 0){
      z_score <- add_NA(z_score, seed_position)
    }
    fwer_critical_value <- hypo$fwer_critical_value
    
    # regions with significant associations
    n_fwer <- length(fwer)
    if(n_fwer == 1){
      significant_regions <- (abs(z_score) > fwer_critical_value)
    }else{
      significant_regions <- matrix(nrow = n_fwer, ncol = length(z_score))
      for(j in 1:n_fwer){
        significant_regions[j,] <- (abs(z_score) > fwer_critical_value[j])
      }
    }
  }
  
  if(length(seed_position) != 0){
    adj_association <- add_NA(adj_association, seed_position)
  }
  
  if(test){
    return(list(
      est = est,
      adj_association = adj_association,
      z_score = z_score,
      fwer_critical_value = fwer_critical_value,
      significant_regions = significant_regions,
      density_log = density_log,
      weight_diagnostics = weight_diagnostics)
    )
  }
  
  return(list(
    est = est,
    adj_association = adj_association,
    density_log = density_log,
    weight_diagnostics = weight_diagnostics)
  )
}

