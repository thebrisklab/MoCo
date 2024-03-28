#' Non-parametric efficient motion-controlled functional connectivity estimators
#' 
#' @details Compute non-parametric efficient motion-controlled functional connectivity estimators (MoCo)
#' 
#' @param A Binary vector of length n \code{number of participants}, denoting diagnosis status.
#' @param X Dataframe or matrix of baseline covariates, typically demographic covariates that would be balanced in a randomized control trial.
#' @param Z Dataframe or matrix of diagnosis-related covariates.
#' @param M Numeric vector of length n representing the motion values.
#' @param Y Matrix of dimension n \times p, where n is the number of participants, and p is the number of functional connectivity of interest. 
#'          Each (i, j) element denotes participant i's functional connectivity between the seed node and node j.
#' @param Delta_M Binary vector of length n indicating data usability, based on the motion value. 
#'                It corresponds to a binary variable indicating whether motion is available and meets inclusion criteria for conventional analyses.
#' @param thresh Value used to threshold M to produce Delta_M. One can specify either Delta_M or thresh.
#' @param Delta_Y Binary vector indicating the quality of T1-weighted image or missingness of rs-fMRI data. 
#' 
#' @param SL_library SuperLearner library for estimating nuisance regressions. Defaults to SL_library if not specified.
#' @param glm_formula All glm formulas default to NULL, indicating SuperLearner will be used for nuisance regressions.
#'                    - \code{gA}: GLM formula for estimating the propensity score.
#'                    - \code{gDM}: GLM formula for estimating the probability P(Delta_M = 1 | A, X).
#'                    - \code{gDY_AX}: GLM formula for estimating the probability P(Delta_Y = 1 | A, X).
#'                    - \code{gDY_AXZ}: GLM formula for estimating the probability P(Delta_Y = 1 | A, X, Z).
#'                    - \code{mu_AMXZ}: GLM formula for estimating the outcome regression E(Y | Delta_Y = 1, A, M, X, Z).
#'                    - \code{eta_AXZ}: GLM formula for estimating E(mu_AMXZ pMXD / pMXZD | A, X, Z, Delta_M = 1).                
#'                    - \code{eta_AXM}: GLM formula for estimating E(mu_AMXZ pMX/pMXZ gDY_AX/gDY_AXZ | A, M, X, Delta_Y = 1).
#'                    - \code{xi_AX}: GLM formula for estimating E(eta_AXZ | A, X).
#'                    - \code{pMX}: GLM formula for estimating p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M = 1), assuming M follows a normal distribution.
#'                    - \code{pMXZ}: GLM formula for estimating p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M = 1), assuming M follows a normal distribution.
#'                    
#' @param HAL_pMX Specifies whether to estimate p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M=1) using the highly adaptive lasso conditional density estimation method. Defaults to \code{TRUE}. 
#' @param HAL_pMXZ Specifies whether to estimate p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M=1) using the highly adaptive lasso conditional density estimation method. Defaults to \code{TRUE}. 
#' 
#' @param HAL_options Additional options for highly adaptive lasso (HAL) method.
#'                   - \code{max_degree}: The highest order of interaction terms for generating basis functions (passed to \code{haldensify}).
#'                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression (passed to \code{haldensify}).
#'                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions (passed to \code{haldensify}).
#' 
#' @param cross_fit Logical indicating whether to develop the estimator based on cross-fitting. Defaults to TRUE.
#' @param test Logical indicating whether to conduct hypothesis testing based on simultaneous confidence band. Defaults to TRUE.
#' @param fwer Family-wise error rate (FWER) to control for multiple hypothesis testing. Defaults to 0.05. Set to NULL if hypo_test is FALSE.
#' @param seed Specifies the value of seed(s) for nuisance regression calculation using super learner. Can be a vector. Defaults to value 1. 
#' 
#' @importFrom SuperLearner SuperLearner
#' @importFrom haldensify haldensify
#' 
#' @return A list with named entries 
#' \describe{
#'   \item{est}{A two times p matrix showing the one-step estimators of the control group and the disease group for each functional connectivity of interest, respectively.}
#'   \item{adj_association}{A p-length vector showing the motion-adjusted association for each functional connectivity of interest, respectively.}
#'   \item{z_score}{A p-length vector of z_score for each functional connectivity.}
#'   \item{significant_regions}{A p-length vector of TRUE or FALSE, indicating significance for each functional connectivity.}
#' }
#' 
#' @export

moco <- function(
  X, Z, A, M, Y, 
  Delta_M, 
  thresh = NULL,
  Delta_Y,
  SL_library = c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.glm.interaction","SL.ranger", "SL.xgboost","SL.mean"),
  glm_formula = list(gA = NULL, 
                     gDM = NULL,
                     gDY_AX = NULL,
                     gDY_AXZ = NULL,
                     mu_AMXZ = NULL,
                     eta_AXZ = NULL,
                     xi_AX = NULL,
                     eta_AXM = NULL,
                     pMX = NULL,
                     pMXZ = NULL),
  HAL_pMX = TRUE,
  HAL_pMXZ = TRUE,
  HAL_options = list(max_degree = 6, lambda_seq = exp(seq(-1, -10, length = 100)), num_knots = c(1000, 500, 250)),
  cross_fit = TRUE,
  test = TRUE,
  fwer = 0.05, 
  seed = 1, 
  ...
){
  
  # total number of runs for handling Monte Carlo variability
  r <- length(seed)
  
  # create matrix and vector for storing fitting results
  est <- list()
  adj_association <- matrix(nrow = r, ncol = ncol(Y))
  
  # create matrix and vector for storing information from hypothesis testing
  if(test){
    # matrix for storing z_score
    z_score_mat <- matrix(nrow = r, ncol = ncol(Y))
    # vector for storing confidence band
    conf_band <- numeric(r)
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
        glm_formula,
        HAL_pMX = HAL_pMX,
        HAL_pMXZ = HAL_pMXZ,
        HAL_options = HAL_options,
        seed = seed[i], 
        cv_folds = cv_folds
      )
    }else{
      result <- one_step(
        X, Z, A, M, Y, 
        Delta_M, 
        thresh,
        Delta_Y,
        SL_library,
        glm_formula,
        HAL_pMX = HAL_pMX,
        HAL_pMXZ = HAL_pMXZ,
        HAL_options = HAL_options,
        seed = seed[i]
      )
    }
    
    # store estimation results
    est <- c(est, list(result$est))
    adj_association[i,] <- result$adj_association
    
    # conduct hypothesis testing based on simultaneous confidence band
    if(test){
      
      hypo <- hypo_test(result, fwer = fwer)
      
      # store z_score and confidence band
      z_score_mat[i,] <- hypo$z_score
      conf_band[i] <- hypo$conf_band
    }
  }
  
  # merge multi-run estimation results
  est <- Reduce("+", est)/r
  est <- as.data.frame(est, row.names = c('est_A0', 'est_A1'))
  colnames(est) <- NULL
  
  # motion-controlled association
  adj_association <- colMeans(adj_association)
  
  if(test){
    # calculate z_score average across runs
    z_score <- colMeans(z_score_mat)
    # calculate average confidence band
    conf_band_avg <- mean(conf_band)
    
    # regions with significant associations
    significant_regions <- (abs(z_score) > conf_band_avg)
  }

  return(list(
    est = est,
    adj_association = adj_association,
    z_score = z_score,
    significant_regions = significant_regions)
  )
}

