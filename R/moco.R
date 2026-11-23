#' Estimate motion-controlled group means and associations
#' 
#' @details Compute motion-controlled brain phenotype estimators (MoCo): estimate group-specific brain-phenotype means and their contrast using one-step estimation, optional cross-fitting and simulations inference.
#' 
#' @param X A dataframe or matrix containing demographic confounders that would ideally be balanced in a randomized controlled trial.
#' @param Z A dataframe or matrix of covariates representing brain phenotypes.
#' @param A A binary vector of length n (number of participants), serving as a group indicator, such as diagnosis group (1) or control group (0).
#' @param M A numeric vector of length n representing continuous motion values for each participant. For example mean FD in mm.
#'          There is no default: set \code{M = NULL} explicitly when a continuous motion measure is unavailable. \code{Delta_M} is then still required
#'          (\code{thresh} cannot be used), usable data are redefined as \code{Delta_Y = 1} and \code{Delta_M = 1}, and participants who fail either
#'          criterion are treated as missing data, similar to Nebel et al. (2022). No motion densities are fitted, so \code{density_method} and its options are not used.
#' @param Y A matrix of dimension n x p, where n is the number of participants, and p is the number of regions of interest.
#'          If it represents seed-based association measures:
#'            Each (i, j) element denotes participant i's association measure between the seed region and region j.
#'            The column representing the association measure of the seed region with itself should be filled with NA values to indicate its position.
#'          If it represents other types of association measures:
#'            Each (i, j) element denotes participant i's association measure between two brain regions of interest, such as the upper diagonal part of the functional connectivity matrix.
#'            No NA values are allowed in \code{Y} in this case.
#' @param Delta_M A binary (0/1) vector of length n indicating whether motion is available and meets the motion quality control (QC) inclusion criteria.
#'          Set \code{Delta_M} = 1 if the participant passes motion QC; otherwise, set \code{Delta_M} = 0.
#'          Any QC definition can be used, e.g., a threshold on mean framewise displacement (\code{Delta_M = as.numeric(M < 0.2)}),
#'          or criteria based on scrubbing or the number of remaining frames. \code{Delta_M} must be 0 when \code{M} is missing.
#'          Required unless \code{thresh} is given (then \code{Delta_M = NULL}, the default, defines it from \code{thresh}).
#' @param Delta_Y A binary vector indicating the non-missingness and whether the brain image data \code{Y} passes quality control after preprocessing. 
#'          Set \code{Delta_Y = 1} if \code{Y} is usable; otherwise, set \code{Delta_Y = 0}.
#' @param thresh \code{NULL} (default) or one positive number c for a motion threshold, \code{Delta_M = 1(M <= c)}.
#'          When \code{Delta_M = NULL}, \code{Delta_M} is defined as \code{as.numeric(!is.na(M) & M <= thresh)}; when both are given,
#'          \code{Delta_M} must equal that rule. With \code{density_method = "gamlss"} or \code{"glm"}, \code{pMXD} and \code{pMXZD} are then fitted as
#'          untruncated GG or log-normal densities on all participants with observed \code{M} and truncated analytically at c,
#'          p(m | Delta_M = 1, a, x) = f(m | a, x) / F(c | a, x) for m <= c and 0 above c, with z added for \code{pMXZD}
#'          (f and F are the fitted density and CDF; F(c | a, x) is floored at 0.001).
#'          \code{pMX} and \code{pMXZ} are fitted on \code{Delta_Y = 1} as usual. With \code{"hal"}, \code{thresh} only defines \code{Delta_M}.
#' 
#' @param SL_library SuperLearner library for estimating nuisance regressions. 
#'                   Defaults to c("SL.earth","SL.glmnet","SL.mgcv_gam","SL.glm", "SL.glm.interaction", "SL.step","SL.step.interaction","SL.xgboost","SL.ranger","SL.mean") if not specified.
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
#'          \code{"gamlss"} (generalized gamma GAMLSS, controlled by \code{gamlss_options}).
#'          Upper or lower case is accepted. A length-2 vector \code{c(pMX = , pMXZ = )} sets the density pairs separately.
#'          Options for a method that is not used are ignored, with a message.
#' @param HAL_options Options for \code{density_method = "hal"}, passed to \code{haldensify::haldensify()}:
#'                   - \code{max_degree}: The highest order of interaction terms for generating basis functions.
#'                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression.
#'                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions.
#' @param glm_options Right-hand sides of the log-normal GLM motion densities, used only when \code{density_method = "glm"}:
#'                   - \code{pMX}: formula for p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M = 1). Default \code{"."} (all of \code{A} and \code{X}).
#'                   - \code{pMXZ}: formula for p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M = 1). Default \code{"."} (all of \code{A}, \code{X} and \code{Z}).
#'          With \code{thresh}, \code{pMXD} and \code{pMXZD} are truncated log-normal densities (see \code{thresh}).
#' @param gamlss_options Options for \code{density_method = "gamlss"}, which fits generalized gamma (GG) densities with \code{gamlss::gamlss()}
#'          (\code{gamlss} and \code{gamlss.dist} are installed with MoCo). Entries left out take their defaults:
#'                   - \code{mu}: \code{"smooth"} (continuous covariates entered as \code{pb()} P-splines, the others linearly; default), \code{"linear"} (all covariates entered linearly), or a one-sided formula.
#'                   - \code{sigma}: \code{"linear"} (default), \code{"smooth"}, \code{"1"} (constant), \code{"A"} (group only), or a one-sided formula.
#'                   - \code{nu}: \code{"1"} (constant nu, estimated; default) or \code{0} (nu fixed at 0, i.e. the log-normal distribution).
#'                   - \code{smooth}: names of the continuous covariates that get \code{pb()} when \code{mu} or \code{sigma} is \code{"smooth"}.
#'                     \code{NULL} (default) uses every numeric covariate with more than 10 distinct values. A \code{pb()} term on a covariate that is mostly zeros
#'                     (for example a score that is 0 for one group) can fail to converge; name the covariates to smooth in that case.
#'                   - \code{method}: \code{"RS"} (default), \code{"CG"}, or \code{"mixed"}. Mixed calls \code{mixed()} with GAMLSS defaults: one RS iteration followed by up to 20 CG iterations, independent of \code{n.cyc}.
#'                   - \code{n.cyc}: maximum outer iterations for RS or CG (default 500); a positive whole number. It does not set the mixed schedule.
#'                   - \code{c.crit}: positive convergence tolerance (default 0.001).
#'                   - \code{trace}: \code{TRUE} (default) prints the global deviance at each GAMLSS iteration (\code{gamlss.control(trace = )}); \code{FALSE} fits silently.
#'          With \code{thresh}, \code{pMXD} and \code{pMXZD} are truncated GG densities (see \code{thresh}).
#'          Each density is fitted once using the selected algorithm and controls (defaults: RS, 500 cycles, tolerance 0.001). 
#' @param cross_fit Logical indicating whether to develop the estimator based on cross-fitting. Defaults to TRUE.
#' @param cv_folds Number of folds used for cross-fitting when \code{cross_fit = TRUE}. Defaults to 5.
#' @param seed_rng Specifies the value of seed(s) for nuisance regression calculation using super learner. Can be a vector. Defaults to value 1. When several seeds are given, the full analysis is run once per seed; the estimates and the participant-level efficient influence functions are averaged across seeds, the covariance is recomputed from the averaged efficient influence functions, and one simultaneous test is conducted across all outcomes in \code{Y}.
#' @param test Logical indicating whether to conduct hypothesis testing based on simultaneous confidence band. Defaults to TRUE.
#' @param fwer A vector of family-wise error rates (FWER) to control for multiple hypothesis testing. Defaults to c(0.05). Set to NULL if test is FALSE.
#' @param hypo_method Critical-value simulation: "multiplier" (Gaussian multiplier bootstrap, default) or "mvn" (multivariate normal draws with the EIF correlation matrix). See \code{\link{hypo_test}}.
#' @param hypo_n_sim Number of Monte Carlo draws for the simultaneous test. Defaults to 10000.
#' @param hypo_seed Critical-value simulation seed, separate from nuisance seeds.
#' 
#' @import SuperLearner 
#' @import haldensify 
#' @import MASS
#' @references Nebel, M. B., Lidstone, D. E., Wang, L., Benkeser, D., Mostofsky, S. H. and Risk, B. B. (2022).
#'   Accounting for motion in resting-state fMRI: What part of the spectrum are we characterizing in autism spectrum disorder?
#'   \emph{NeuroImage}, 257, 119296. \doi{10.1016/j.neuroimage.2022.119296}
#' 
#' @importFrom stats binomial gaussian dnorm sd cov
#' 
#' @return A list with named entries 
#' \describe{
#'   \item{est}{A two times p matrix showing the one-step estimators of the control group and the disease group for each functional connectivity of interest, respectively.}
#'   \item{adj_association}{A p-length vector showing the motion-adjusted association for each functional connectivity of interest, respectively.}
#'   \item{z_score}{A p-length vector of z_score for each functional connectivity if \code{test} is TRUE.}
#'   \item{fwer_critical_value}{Critical value(s) of the simultaneous confidence band (the (1 - FWER) quantile of the maximum absolute z-statistic) for each specified FWER.}
#'   \item{significant_regions}{A p-length vector of TRUE or FALSE, indicating significance for each functional connectivity if \code{test} is TRUE.}
#' }
#' 
#' @examples
#' data(moco_example)
#' d <- moco_example
#'
#' # Fast example (a few seconds): one outcome (Frontoparietal) plus the
#' # seed-to-self column, a small Super Learner library, log-normal GLM motion
#' # densities, and 10,000 draws for the critical value
#' fit <- moco(X = d$X, Z = d$Z, A = d$A, M = d$M, Y = d$Y[, 6:7],
#'             Delta_M = d$Delta_M, Delta_Y = d$Delta_Y,
#'             SL_library = c("SL.mean", "SL.glm"),
#'             density_method = "glm", cv_folds = 5, seed_rng = 1,
#'             hypo_n_sim = 10000)
#' fit$adj_association
#' fit$significant_regions
#'
#' \donttest{
#' # The same analysis with GAMLSS motion densities (about 20 s): smooth terms
#' # for age (X2) and IQ (Z2)
#' fit_gamlss <- moco(X = d$X, Z = d$Z, A = d$A, M = d$M, Y = d$Y[, 6:7],
#'                    Delta_M = d$Delta_M, Delta_Y = d$Delta_Y,
#'                    SL_library = c("SL.mean", "SL.glm"),
#'                    density_method = "gamlss",
#'                    gamlss_options = list(smooth = c("X2", "Z2")),
#'                    cv_folds = 5, seed_rng = 1)
#' fit_gamlss$adj_association
#' }
#'
#' # The full GLM, GAMLSS and HAL analysis of the README tutorial (about 1 hour
#' # in total) is in an installed script:
#' # file.show(system.file("examples", "full_analysis.R", package = "MoCo"))
#'
#' @export

moco <- function(
  X, Z, A, M, Y, 
  Delta_M = NULL, 
  Delta_Y,
  thresh = NULL,
  SL_library = c("SL.earth","SL.glmnet","SL.mgcv_gam","SL.glm", "SL.glm.interaction", "SL.step","SL.step.interaction","SL.xgboost","SL.ranger","SL.mean"),
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
  glm_formula = list(
    gA = NULL, 
    gDM = NULL,
    gDY_AX = NULL,
    gDY_AXZ = NULL,
    mu_AMXZ = NULL,
    eta_AXZ = NULL,
    eta_AXM = NULL,
    xi_AX = NULL
  ),
  density_method = "hal",
  HAL_options = list(
    max_degree = 3, 
    lambda_seq = exp(seq(-1, -10, length = 100)), 
    num_knots = c(1000, 500, 250)
    ),
  glm_options = list(
    pMX = ".", 
    pMXZ = "."
  ),
  gamlss_options = list(
  mu = "smooth",
  sigma = "linear",
  nu = "1",
  smooth = NULL,
  method = "RS",
  n.cyc = 500,
  c.crit = 0.001,
  trace = TRUE
  ),
  cross_fit = TRUE,
  cv_folds = 5,
  seed_rng = 1, 
  test = TRUE,
  fwer = 0.05, 
  hypo_method = "multiplier",
  hypo_n_sim = 10000L,
  hypo_seed = 1
){
  
  # Leave the caller's random-number state as it was before the call.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  
  # which density options the user supplied 
  HAL_options_given <- !missing(HAL_options)
  glm_options_given <- !missing(glm_options)
  gamlss_options_given <- !missing(gamlss_options)
  
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
  
  # thresh: Delta_M = 1(M <= thresh); the GAMLSS and GLM densities pMXD/pMXZD are then truncated at thresh
  if(!is.null(thresh)){
    if(!is.numeric(thresh) || length(thresh) != 1 || !is.finite(thresh) || thresh <= 0){
      stop("thresh must be one positive number, e.g. thresh = 0.2.")
    }
    if(is.null(M) || length(M) != n){
      stop("thresh requires M, a motion vector of length n = nrow(Y).")
    }
    Delta_M_thresh <- as.numeric(!is.na(M) & M <= thresh)
    if(is.null(Delta_M)){
      Delta_M <- Delta_M_thresh
    }else if(length(Delta_M) != n || any(is.na(Delta_M)) || any(Delta_M != Delta_M_thresh)){
      stop("Delta_M is not equal to 1(M <= thresh). When both are given, Delta_M must be the motion threshold rule; ",
           "set Delta_M = NULL to define it from thresh, or thresh = NULL to use Delta_M without truncation.")
    }
  }
  
  # check Delta_M: required for every definition of motion quality control
  if(is.null(Delta_M)){
    stop("Delta_M or thresh is required. Please provide a 0/1 vector indicating whether each participant passes motion quality control, ",
         "or thresh = c to define Delta_M = 1(M <= c).")
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
  
  # motion density method
  dens_method <- resolve_density_method(density_method)
  if(!is.null(M)){
    if(HAL_options_given && !any(dens_method == "hal")){
      message("HAL_options is only used by density_method = \"hal\"; it is ignored here.")
    }
    if(glm_options_given && !any(dens_method == "glm")){
      message("glm_options is only used by density_method = \"glm\"; it is ignored here.")
    }
    if(gamlss_options_given && !is.null(gamlss_options) && !any(dens_method == "gamlss")){
      message("gamlss_options is only used by density_method = \"gamlss\"; it is ignored here.")
    }
  }
  uses_gamlss <- !is.null(M) && any(dens_method == "gamlss")
  if(uses_gamlss){
    gg_check_packages()
    gg_check_options(gamlss_options)
  }
  # Delta_M that looks like a motion threshold, but no thresh given (GAMLSS and GLM densities)
  if(!is.null(M) && any(dens_method %in% c("gamlss", "glm"))){
    obs <- !is.na(M)
    if(is.null(thresh) && any(obs & Delta_M == 1) && any(obs & Delta_M == 0) &&
       max(M[obs & Delta_M == 1]) < min(M[obs & Delta_M == 0])){
      message("Delta_M looks like a motion threshold: all M with Delta_M = 1 are below all M with Delta_M = 0. ",
              "If Delta_M = 1(M <= c), supply thresh = c to use the truncated GAMLSS or GLM densities.")
    }
  }
  
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
  
  for(i in 1:r){
    
    # whether to use cross-fit or not
    if(cross_fit){
      result <- one_step_cross(
        X = X, Z = Z, A = A, M = M, Y = Y, 
        Delta_M = Delta_M, 
        Delta_Y = Delta_Y,
        thresh = thresh,
        SL_library = SL_library,
        SL_library_customize = SL_library_customize,
        glm_formula = glm_formula,
        HAL_options = HAL_options,
        glm_options = glm_options,
        density_method = dens_method,
        gamlss_options = gamlss_options,
        seed = seed_rng[i], 
        cv_folds = cv_folds
      )
    }else{
      result <- one_step(
        X = X, Z = Z, A = A, M = M, Y = Y, 
        Delta_M = Delta_M, 
        Delta_Y = Delta_Y,
        thresh = thresh,
        SL_library = SL_library,
        SL_library_customize = SL_library_customize,
        glm_formula = glm_formula,
        HAL_options = HAL_options,
        glm_options = glm_options,
        density_method = dens_method,
        gamlss_options = gamlss_options,
        seed = seed_rng[i]
      )
    }
    
    # store estimation results
    est <- c(est, list(result$est))
    adj_association[i,] <- result$adj_association
    # accumulate participant-level efficient influence functions for hypothesis testing
    if(test){
      if(is.null(eif_sum)){
        eif_sum <- result$eif_mat
      }else{
        eif_sum <- Map(`+`, eif_sum, result$eif_mat)
      }
    }
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
      fwer = fwer, method = hypo_method,
      n_sim = hypo_n_sim, seed = hypo_seed
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
      significant_regions = significant_regions)
    )
  }
  
  return(list(
    est = est,
    adj_association = adj_association)
  )
}

