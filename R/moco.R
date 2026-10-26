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
#' @param Delta_M A binary vector of length n indicating whether motion is available and meets inclusion criteria. 
#'          If motion meets inclusion criteria for analysis, set \code{Delta_M} = 1; otherwise, set \code{Delta_M} = 0.             
#' @param thresh A numeric value used to threshold M to produce \code{Delta_M}. One can specify either \code{Delta_M} or thresh.
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
#'                    - \code{pMX}: GLM formula for estimating p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M = 1), assuming M follows a log normal distribution.
#'                    - \code{pMXZ}: GLM formula for estimating p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M = 1), assuming M follows a log normal distribution.
#'                    
#' @param HAL_pMX Specifies whether to estimate p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M=1) using the highly adaptive lasso conditional density estimation method. 
#' Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMX} option in \code{glm_formula}, such as \code{pMX = "."}.
#' @param HAL_pMXZ Specifies whether to estimate p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M=1) using the highly adaptive lasso conditional density estimation method. 
#' Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMXZ} option in \code{glm_formula}, such as \code{pMXZ = "."}.
#'  
#' 
#' @param HAL_options Additional options for highly adaptive lasso (HAL) method. These will be passed to the haldensify function in the haldensify package.
#'                   - \code{max_degree}: The highest order of interaction terms for generating basis functions.
#'                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression.
#'                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions.
#' 
#' @param gamlss_options Optional GAMLSS conditional-density specification (requires the gamlss and gamlss.dist packages). Defaults to \code{NULL} (HAL/GLM as before).
#'                   When supplied, GAMLSS replaces HAL/GLM for all four motion densities: pMX = p(m | a, x) and pMXZ = p(m | a, x, z) are fitted on \code{Delta_Y = 1};
#'                   pMXD = p(m | Delta_M = 1, a, x) and pMXZD = p(m | Delta_M = 1, a, x, z) are fitted on \code{Delta_M = 1}. pMXD (at A = 0) and pMXZD (at the observed A) are evaluated at every observed M,
#'                   with no zeroing where \code{Delta_M = 0}; pMXZ is evaluated at the observed A, A = 0 and A = 1 on \code{Delta_Y = 1}. Each density is fitted once per training fold. A list with (all optional):
#'                   - \code{family}: a gamlss.dist family name (e.g. "LOGNO", "GG", "BCPEo"), family function or family object, or a list by role (\code{pMX}, \code{pMXD}, \code{pMXZ}, \code{pMXZD}) with \code{default}. It must be continuous on (0, Inf); default "GG". All observed M must be > 0.
#'                   - \code{formulas}: list of \code{mu}, \code{sigma}, \code{nu}, \code{tau}, each a formula, a string, or a keyword: "linear" (A, X and, for pMXZ/pMXZD, Z entering linearly), "pb" (the same with \code{pb()} on numeric columns with at least \code{continuous_min_unique} = 10 distinct values) or "1" (intercept only). A "." expands to all available columns. Per-role lists are allowed, e.g. \code{list(mu = "pb", pMXD = list(sigma = "1"))}. Defaults: mu "linear", others "1".
#'                   - \code{method}: gamlss fitting algorithm, e.g. "RS()" (default), "CG()", "mixed(1, 50)" (string or call).
#'                   - \code{control}: list of \code{gamlss.control()} arguments (e.g. \code{n.cyc}, \code{c.crit}, \code{trace}); \code{n.cyc} defaults to 200 and \code{trace} to FALSE. \code{i.control}: list of \code{glim.control()} arguments.
#'                   - \code{weights}, \code{mu.start}, \code{sigma.start}, \code{nu.start}, \code{tau.start}: scalars or length-n vectors (subset to each fold's fit rows); \code{mu.fix}, \code{sigma.fix}, \code{nu.fix}, \code{tau.fix}: logical.
#'                   - \code{args}: named list of any further arguments forwarded to \code{gamlss::gamlss()}.
#'                   - \code{start}: "auto" (default; Box-Cox families BCCG(o)/BCPE(o)/BCT(o) start from the LOGNO fit) or "default". \code{retry}: TRUE (default; on failure or non-convergence refit from the sigma ~ 1 fit, then with CG(), then with RS() if another method was requested) or FALSE.
#'                   - \code{clamp}: TRUE (default; continuous covariates of prediction rows are clamped to that density's training range) or FALSE. \code{on_nonconvergence}: "stop" (default) or "warn".
#'                   - \code{pMX}, \code{pMXD}, \code{pMXZ}, \code{pMXZD}: per-role override lists with any of the entries above, e.g. \code{pMXD = list(family = "BCTo", method = "CG()", control = list(n.cyc = 400))}.
#'                   Per-fold, per-role fit summaries (family, formulas, deviance, AIC, convergence, method, clamping) are returned as \code{gamlss_info}.
#' @param cross_fit Logical indicating whether to develop the estimator based on cross-fitting. Defaults to TRUE.
#' @param cv_folds Number of folds used for cross-fitting when \code{cross_fit = TRUE}. Defaults to 5.
#' @param ... Currently unused; additional arguments are ignored.
#' @param test Logical indicating whether to conduct hypothesis testing based on simultaneous confidence band. Defaults to TRUE.
#' @param fwer A vector of family-wise error rates (FWER) to control for multiple hypothesis testing. Defaults to c(0.05). Set to NULL if test is FALSE.
#' @param seed_rng Specifies the value of seed(s) for nuisance regression calculation using super learner. Can be a vector. Defaults to value 1. 
#' 
#' @import SuperLearner 
#' @import haldensify 
#' @import MASS
#' 
#' @return A list with named entries 
#' \describe{
#'   \item{est}{A two times p matrix showing the one-step estimators of the control group and the disease group for each functional connectivity of interest, respectively.}
#'   \item{adj_association}{A p-length vector showing the motion-adjusted association for each functional connectivity of interest, respectively.}
#'   \item{z_score}{A p-length vector of z_score for each functional connectivity if \code{test} is TRUE.}
#'   \item{fwer_critical_value}{Simultaneous confidence band based on specified FWER.}
#'   \item{gamlss_info}{Only when \code{gamlss_options} is supplied: a data frame with one row per seed, fold and density role summarising the GAMLSS fits.}
#'   \item{significant_regions}{A p-length vector of TRUE or FALSE, indicating significance for each functional connectivity if \code{test} is TRUE.}
#' }
#' 
#' @examples
#' data(data)
#' # Fast illustration: log-normal GLM motion densities and a small library
#' fit <- moco(
#'   X = data$X, Z = data$Z, A = data$A, M = data$M, Y = data$Y,
#'   Delta_M = data$Delta_M, Delta_Y = data$Delta_Y,
#'   SL_library = c("SL.mean", "SL.glm"),
#'   glm_formula = list(pMX = ".", pMXZ = "."),
#'   HAL_pMX = FALSE, HAL_pMXZ = FALSE,
#'   cross_fit = FALSE, test = TRUE, fwer = 0.05
#' )
#' fit$adj_association
#' fit$significant_regions
#'
#' \donttest{
#' # GAMLSS motion densities: generalized gamma by default, log-normal for the
#' # two densities fitted on Delta_M = 1
#' fit_gamlss <- moco(
#'   X = data$X, Z = data$Z, A = data$A, M = data$M, Y = data$Y,
#'   Delta_M = data$Delta_M, Delta_Y = data$Delta_Y,
#'   SL_library = c("SL.mean", "SL.glm"),
#'   gamlss_options = list(
#'     family = list(default = "GG", pMXD = "LOGNO", pMXZD = "LOGNO")
#'   ),
#'   cross_fit = FALSE, test = TRUE, fwer = 0.05
#' )
#' fit_gamlss$adj_association
#' }
#'
#' @export

moco <- function(
  X, Z, A, M, Y, 
  Delta_M, 
  thresh = NULL,
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
                     xi_AX = NULL,
                     pMX = NULL,
                     pMXZ = NULL),
  HAL_pMX = TRUE,
  HAL_pMXZ = TRUE,
  HAL_options = list(max_degree = 3, lambda_seq = exp(seq(-1, -10, length = 100)), num_knots = c(1000, 500, 250)),
  cross_fit = TRUE,
  cv_folds = 5,
  test = TRUE,
  fwer = 0.05, 
  seed_rng = 1, 
  gamlss_options = NULL,
  ...
){
  
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
  
  # seed position
  seed_position <- which(apply(Y, 2, function(x){
    sum(is.na(x))
  }) == nrow(Y))
  
  if(length(seed_position) != 0){
    # remove NA column for subsequent calculation
    Y <- Y[, -seed_position]
  }
  
  # create matrix and vector for storing fitting results
  est <- list()
  gamlss_info <- list()
  adj_association <- matrix(nrow = r, ncol = ncol(Y))
  
  # create matrix and vector for storing information from hypothesis testing
  if(test){
    # matrix for storing z_score
    z_score_mat <- matrix(nrow = r, ncol = ncol(Y))
    # number of fwer under consideration
    n_fwer <- length(fwer)
    # vector for storing confidence band
    fwer_critical_value <- matrix(nrow = r, ncol = n_fwer)
  }
  
  for(i in 1:r){
    
    # whether to use cross-fit or not
    if(cross_fit){
      result <- one_step_cross(
        X, Z, A, M, Y, 
        Delta_M, 
        thresh,
        Delta_Y,
        SL_library,
        SL_library_customize,
        glm_formula,
        HAL_pMX,
        HAL_pMXZ,
        HAL_options,
        seed = seed_rng[i], 
        cv_folds,
        gamlss_options = gamlss_options
      )
    }else{
      result <- one_step(
        X, Z, A, M, Y, 
        Delta_M, 
        thresh,
        Delta_Y,
        SL_library,
        SL_library_customize,
        glm_formula,
        HAL_pMX,
        HAL_pMXZ,
        HAL_options,
        seed = seed_rng[i],
        gamlss_options = gamlss_options
      )
    }
    
    # store estimation results
    est <- c(est, list(result$est))
    adj_association[i,] <- result$adj_association
    if(!is.null(result$gamlss_info)) gamlss_info[[i]] <- cbind(seed = seed_rng[i], result$gamlss_info)
    
    # conduct hypothesis testing based on simultaneous confidence band
    if(test){
      
      hypo <- hypo_test(result, fwer = fwer)
      
      # store z_score and confidence band
      z_score_mat[i,] <- hypo$z_score
      fwer_critical_value[i,] <- hypo$fwer_critical_value
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
  if(length(seed_position) != 0){
    adj_association <- add_NA(adj_association, seed_position)
  }
  
  if(test){
    # calculate z_score average across runs
    z_score <- colMeans(z_score_mat)
    if(length(seed_position) != 0){
      z_score <- add_NA(z_score, seed_position)
    }
    # calculate average confidence band
    fwer_critical_value_avg <- colMeans(fwer_critical_value)
    
    # regions with significant associations
    if(n_fwer == 1){
      significant_regions <- (abs(z_score) > fwer_critical_value_avg)
    }else{
      significant_regions <- matrix(nrow = n_fwer, ncol = length(z_score))
      for(j in 1:n_fwer){
        significant_regions[j,] <- (abs(z_score) > fwer_critical_value_avg[j])
      }
    }
    
    out <- list(
      est = est,
      adj_association = adj_association,
      z_score = z_score,
      fwer_critical_value = fwer_critical_value_avg,
      significant_regions = significant_regions)
    if(length(gamlss_info)) out$gamlss_info <- do.call(rbind, gamlss_info)
    return(out)
  }
  
  out <- list(
    est = est,
    adj_association = adj_association)
  if(length(gamlss_info)) out$gamlss_info <- do.call(rbind, gamlss_info)
  return(out)
}

