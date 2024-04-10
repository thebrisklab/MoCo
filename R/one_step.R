#' Non-parametric efficient motion-controlled functional connectivity estimators
#' 
#' @details Compute non-parametric efficient motion-controlled functional connectivity estimators based on the proposed method.
#' 
#' @param A Binary vector of length n \code{number of participants}, denoting diagnosis status
#' @param X Dataframe or matrix of baseline covariates, typically demographic covariates that would be balanced in a randomized control trial
#' @param Z Dataframe or matrix of diagnosis-related covariates
#' @param M Numeric vector of length n representing the motion values.
#' @param Y Matrix of dimension n \times p, where p is the number of functional connectivity of interest. 
#'          Each (i, j) element denotes participant i's functional connectivity between the seed node and node j.
#' @param Delta_M Binary vector of length n indicating data usability, based on the motion value. 
#'                It corresponds to a binary variable indicating whether motion is available and meets inclusion criteria for conventional analyses.
#' @param thresh Value used to threshold M to produce Delta_M. One can specify either Delta_M or thresh.
#' @param Delta_Y Binary vector indicating the quality of T1-weighted image or missingness of rs-fMRI data. 
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
#' @param HAL_pMX Specifies whether to estimate p(m | a, x) and p(m | a, x, Delta_M=1) using the highly adaptive lasso conditional density estimation method. Defaults to \code{TRUE}. 
#' @param HAL_pMXZ Specifies whether to estimate p(m | a, x, z) and p(m | a, x, z, Delta_M=1) using the highly adaptive lasso conditional density estimation method. Defaults to \code{TRUE}. 
#' 
#' @param HAL_options Additional options for highly adaptive lasso (HAL) method.
#'                   - \code{max_degree}: The highest order of interaction terms for generating basis functions (passed to \code{haldensify}).
#'                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression (passed to \code{haldensify}).
#'                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions (passed to \code{haldensify}).
#' 
#' @importFrom SuperLearner SuperLearner
#' @importFrom haldensify haldensify
#' 
#' @return A list with named entries 
#' \describe{
#'   \item{est}{A two times p matrix showing the one-step estimators of the control group and the disease group for each functional connectivity of interest, respectively.}
#'   \item{adj_association}{A p-length vector showing the motion-controlled association for each functional connectivity of interest, respectively.}
#'   \item{eif_mat}{A p-length list of the estimated EIF evaluated on the observations of the control and disease group, respectively.}
#'   \item{cov_mat}{A p-length list of the estimated covariance matrix of the one-step estimators.}
#' }

one_step <- function(
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
    glm_formula = list(
      gA = NULL, 
      gDM = NULL,
      gDY_AX = NULL,
      gDY_AXZ = NULL,
      mu_AMXZ = NULL,
      eta_AXZ = NULL,
      eta_AXM = NULL,
      xi_AX = NULL,
      pMX = NULL,
      pMXZ = NULL
    ),
    HAL_pMX = TRUE,
    HAL_pMXZ = TRUE,
    HAL_options = list(
      max_degree = 6, 
      lambda_seq = exp(seq(-1, -10, length = 100)), 
      num_knots = c(1000, 500, 250)
    ),
    seed = 1, 
    ...
){
  # number of participants
  n <- nrow(Y)
  # number of edges
  p <- ncol(Y)
  
  if(!is.null(thresh)){
    # if thresh is not null, collapse M into dummy variables based on truncated level thresh
    Delta_M <- as.numeric(M < thresh)
    Delta_M[is.na(M)] <- 0
  }
  
  # change X and Z to a dataframe if is a vector
  if(is.null(dim(X))){X <- data.frame(X = X)}
  if(is.null(dim(Z))){Z <- data.frame(Z = Z)}
  
  # SL options
  if(is.null(SL_library)){
    SL_gA <- SL_library_customize$gA
    SL_gDM <- SL_library_customize$gDM
    SL_gDY_AX <- SL_library_customize$gDY_AX
    SL_gDY_AXZ <- SL_library_customize$gDY_AXZ
    SL_mu_AMXZ <- SL_library_customize$mu_AMXZ
    SL_eta_AXZ <- SL_library_customize$eta_AXZ
    SL_eta_AXM <- SL_library_customize$eta_AXM
    SL_xi_AX <- SL_library_customize$xi_AX
  }else{
    SL_gA <- SL_gDM <- SL_gDY_AX <- SL_gDY_AXZ <- SL_mu_AMXZ <- SL_eta_AXZ <- SL_eta_AXM <- SL_xi_AX <- SL_library
  }

  # fit regression for propensity score 
  # disease status, conditional on baseline covariates
  if(!is.null(glm_formula$gA)){
    gA_fit <- stats::glm(paste0("A ~ ", glm_formula$gA), family = binomial(), 
                         data = data.frame(A = A, X))
    gAn_1 <- stats::predict(gA_fit, type = "response")
  }else{
    set.seed(seed)
    if(ncol(X) == 1){
      SL_gA <- SL_gA[SL_gA != "SL.glmnet"]
    }
    gA_fit <- SuperLearner::SuperLearner(Y = A, X = X,
                                         family = binomial(), 
                                         SL.library = SL_gA,
                                         method = tmp_method.CC_nloglik(),
                                         control = list(saveCVFitLibrary = TRUE))
    gAn_1 <- stats::predict(gA_fit, type = "response", newdata = X)[[1]]
  }
  
  # fit binary mediator regression
  if(!is.null(glm_formula$gDM)){
    gDM_fit <- stats::glm(paste0("Delta_M ~ ", glm_formula$gDM), family = binomial(), 
                          data = data.frame(Delta_M = Delta_M, A, X))
    gDMn_1_A0 <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, X))
  }else{
    set.seed(seed)
    gDM_fit <- SuperLearner::SuperLearner(Y = Delta_M, X = data.frame(A, X),
                                          family = binomial(), 
                                          SL.library = SL_gDM,
                                          method = tmp_method.CC_nloglik(),
                                          control = list(saveCVFitLibrary = TRUE))
    gDMn_1_A0 <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, X))[[1]]
  }
  
  # fit missingness indicator regression
  if(sum(Delta_Y) == n){
    gDYn_1_AX <- rep(1, n)
    gDYn_1_AXZ <- gDYn_1_A0XZ <- gDYn_1_A1XZ <- rep(1, n)
  }else{
    # probability of non-missing conditioning on disease status A, baseline covariates X and diagnosis-related covariates Z
    if(!is.null(glm_formula$gDY_AX)){
      gDY_AX_fit <- stats::glm(paste0("Delta_Y ~ ", glm_formula$gDY_AX), family = binomial(), 
                               data = data.frame(Delta_Y = Delta_Y, A, X))
      gDYn_1_AX <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(A, X))
    }else{
      set.seed(seed)
      gDY_AX_fit <- SuperLearner::SuperLearner(Y = Delta_Y, X = data.frame(A, X),
                                               family = binomial(), 
                                               SL.library = SL_gDY_AX,
                                               method = tmp_method.CC_nloglik(),
                                               control = list(saveCVFitLibrary = TRUE))
      gDYn_1_AX <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(A, X))[[1]]
    }
  
    # probability of non-missing conditioning on disease status A, baseline covariates X and diagnosis-related covariates Z
    if(!is.null(glm_formula$gDY_AXZ)){
      gDY_AXZ_fit <- stats::glm(paste0("Delta_Y ~ ", glm_formula$gDY_AXZ), family = binomial(), 
                                data = data.frame(Delta_Y = Delta_Y, A, X, Z))
      gDYn_1_AXZ <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A, X, Z))
      gDYn_1_A0XZ <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, X, Z))
      gDYn_1_A1XZ <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, X, Z))
    }else{
      set.seed(seed)
      gDY_AXZ_fit <- SuperLearner::SuperLearner(Y = Delta_Y, X = data.frame(A, X, Z),
                                                family = binomial(), 
                                                SL.library = SL_gDY_AXZ,
                                                method = tmp_method.CC_nloglik(),
                                                control = list(saveCVFitLibrary = TRUE))
      gDYn_1_AXZ <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A, X, Z))[[1]]
      gDYn_1_A0XZ <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, X, Z))[[1]]
      gDYn_1_A1XZ <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, X, Z))[[1]]
    }
  }
  
  if(HAL_pMX){
    # density estimation for mediator, conditioning on disease status and baseline covariates p(m|a,x)
    # estimate p(m|a,x) use highly adaptive lasso conditional density estimation method
    # use default n_bins range, use cv to choose number of bins
    pMX_fit <- haldensify::haldensify(
      A = M[Delta_Y == 1],
      W = data.frame(A, X)[Delta_Y == 1,],
      max_degree = HAL_options$max_degree,
      lambda_seq = HAL_options$lambda_seq,
      num_knots = HAL_options$num_knots
    ) 
    pMXn_A <- rep(NA, n)
    pMXn_A[Delta_Y == 1] <- stats::predict(pMX_fit, new_A = M[Delta_Y == 1], new_W = data.frame(A, X)[Delta_Y == 1,])
    
    # density estimation for mediator, conditioning on diease status and baseline covariates p(m|0,x,Delta_M=1)
    pMXD_fit <- haldensify::haldensify(
      A = M[Delta_M == 1],
      W = data.frame(A, X)[Delta_M == 1,],  
      max_degree = HAL_options$max_degree,
      lambda_seq = HAL_options$lambda_seq,
      num_knots = HAL_options$num_knots
    ) 
    pMXDn_A0 <- rep(NA, n)
    pMXDn_A0[which(!is.na(M))] <- stats::predict(pMXD_fit, new_A = M[which(!is.na(M))], new_W = data.frame(A = 0, X)[which(!is.na(M)),], trim_min = 0)
  }else{
    pMX_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M), A, X)[Delta_Y == 1,])
    pMXn_A <- rep(NA, n)
    pMXn_A <- (1/M)[Delta_Y == 1] * dnorm(log(M)[Delta_Y == 1], mean = stats::predict(pMX_fit, newdata = data.frame(A, X)[Delta_Y == 1,]), sd = sd(pMX_fit$residuals)) 
    
    pMXD_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M), A, X)[Delta_M == 1,])
    pMXDn_A0 <- rep(NA, n)
    pMXDn_A0[Delta_M==0] <- 0
    pMXDn_A0[Delta_M==1] <- (1/M)[Delta_M == 1] * dnorm(log(M)[Delta_M == 1], mean = stats::predict(pMXD_fit, newdata = data.frame(A=0,X)[Delta_M == 1,]), sd = sd(pMXD_fit$residuals)) 
  }
  
  if(HAL_pMXZ){
    # density estimation for mediator, conditioning on disease status, binary mediator, 
    # baseline covariates, and mediator-outcome confounder p(m|0,x,z)
    # use default n_bins range, use cv to choose number of bins
    pMXZ_fit <- haldensify::haldensify(
      A = M[Delta_Y == 1],
      W = data.frame(A, X, Z)[Delta_Y == 1,], 
      max_degree = HAL_options$max_degree,
      lambda_seq = HAL_options$lambda_seq,
      num_knots = HAL_options$num_knots
    )
    pMXZn_A0 <- pMXZn_A1 <- pMXZn_A <- rep(NA, n)
    pMXZn_A[Delta_Y == 1] <- stats::predict(pMXZ_fit, new_A = M[Delta_Y == 1], new_W = data.frame(A, X, Z)[Delta_Y == 1,])
    pMXZn_A0[Delta_Y == 1] <- stats::predict(pMXZ_fit, new_A = M[Delta_Y == 1], new_W = data.frame(A = 0, X, Z)[Delta_Y == 1,])
    pMXZn_A1[Delta_Y == 1] <- stats::predict(pMXZ_fit, new_A = M[Delta_Y == 1], new_W = data.frame(A = 1, X, Z)[Delta_Y == 1,])
    
    # density estimation for mediator, conditioning on disease status, binary mediator, 
    # baseline covariates, and mediator-outcome confounder p(m|0,x,z,Delta_M=1)
    pMXZD_fit <- haldensify::haldensify(
      A = M[Delta_M == 1],
      W = data.frame(A, X, Z)[Delta_M == 1,], 
      max_degree = HAL_options$max_degree,
      lambda_seq = HAL_options$lambda_seq,
      num_knots = HAL_options$num_knots
    ) 
    # predict and set probability for M value out of the support to be 0
    pMXZDn_A <- rep(NA, n)
    pMXZDn_A[which(!is.na(M))] <- stats::predict(pMXZD_fit, new_A = M[which(!is.na(M))], new_W = data.frame(A, X, Z)[which(!is.na(M)),], trim_min = 0)
  }else{
    pMXZ_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMXZ), family = gaussian(), data = data.frame(log_M = log(M), A, X, Z)[Delta_Y == 1,])
    pMXZn_A0 <- pMXZn_A1 <- pMXZn_A <- rep(NA, n)
    pMXZn_A[Delta_Y == 1] <- (1/M)[Delta_Y == 1] * dnorm(log(M)[Delta_Y == 1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A, X, Z)[Delta_Y == 1,]), sd = sd(pMXZ_fit$residuals))    
    pMXZn_A0[Delta_Y == 1] <- (1/M)[Delta_Y == 1] * dnorm(log(M)[Delta_Y == 1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 0, X, Z)[Delta_Y == 1,]), sd = sd(pMXZ_fit$residuals))
    pMXZn_A1[Delta_Y == 1] <- (1/M)[Delta_Y == 1] * dnorm(log(M)[Delta_Y == 1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 1, X, Z)[Delta_Y == 1,]), sd = sd(pMXZ_fit$residuals))
    
    pMXZD_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMXZ), family = gaussian(), data = data.frame(log_M = log(M), A, X, Z)[Delta_M == 1,])
    pMXZDn_A <- rep(NA, n)
    pMXZDn_A[Delta_M==0] <- 0
    pMXZDn_A[Delta_M==1] <- (1/M)[Delta_M==1] * dnorm(log(M)[Delta_M==1], mean = stats::predict(pMXZD_fit, newdata = data.frame(A, X, Z)[Delta_M==1,]), sd = sd(pMXZD_fit$residuals))
  }
  
  # define parameters for storing results
  est <- matrix(nrow = 2, ncol = p)
  eif_mat <- list(p)
  cov_mat <- list(p)
  for(j in 1:p){
    # fit outcome regression
    mu_MXZn_A0 <- mu_MXZn_A1 <- mu_MXZn_A <- numeric(n)
    if(!is.null(glm_formula$mu_AMXZ)){
      mu_AMXZ_fit <- stats::glm(paste0("Y ~ ", glm_formula$mu_AMXZ), family = gaussian(), 
                                data = data.frame(Y = Y[,j], A, M, X, Z)[Delta_Y==1,])
      mu_MXZn_A0 <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, M, X, Z))
      mu_MXZn_A1 <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, M, X, Z))
      mu_MXZn_A <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A, M, X, Z))
    }else{
      set.seed(seed)
      mu_AMXZ_fit <- SuperLearner::SuperLearner(Y = Y[Delta_Y==1,j], X = data.frame(A, M, X, Z)[Delta_Y==1,],
                                                family = gaussian(), 
                                                SL.library = SL_mu_AMXZ,
                                                method = tmp_method.CC_LS(),
                                                control = list(saveCVFitLibrary = TRUE))
      mu_MXZn_A0 <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, M, X, Z))[[1]]
      mu_MXZn_A1 <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, M, X, Z))[[1]]
      mu_MXZn_A <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A, M, X, Z))[[1]]
    }
    
    # fit additional pseudo-outcome regression for mu_MXZ*pMXD/pMXZD, conditioning on disease status, binary mediator and baseline covariates
    mu_pseudo_A <- mu_MXZn_A*pMXDn_A0/pMXZDn_A
    if(!is.null(glm_formula$eta_AXZ)){
      eta_AXZ_fit <- stats::glm(paste0("mu_pseudo_A ~ ", glm_formula$eta_AXZ), family = gaussian(), 
                                data = data.frame(mu_pseudo_A = mu_pseudo_A, A, X, Z)[Delta_M==1,])
      eta_AXZn_A <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A, X, Z))
      eta_AXZn_A0 <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, X, Z))
      eta_AXZn_A1 <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, X, Z))
    }else{
      set.seed(seed)
      eta_AXZ_fit <- SuperLearner::SuperLearner(Y = mu_pseudo_A[Delta_M==1], X = data.frame(A, X, Z)[Delta_M==1,],
                                                family = gaussian(),
                                                SL.library = SL_eta_AXZ,
                                                method = tmp_method.CC_LS(),
                                                control = list(saveCVFitLibrary = TRUE))
      eta_AXZn_A <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A, X, Z))[[1]]
      eta_AXZn_A0 <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, X, Z))[[1]]
      eta_AXZn_A1 <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, X, Z))[[1]]
    }
    
    # fit pseudo-outcome regression for Qd, conditioning on disease status, mediator and baseline covariates
    mu_pseudo_A_star <- mu_MXZn_A*pMXn_A/pMXZn_A*gDYn_1_AX/gDYn_1_AXZ
    # fit the regression
    if(!is.null(glm_formula$eta_AXM)){
      eta_AXM_fit <- stats::glm(paste0("mu_pseudo_A_star ~ ", glm_formula$eta_AXM), family = gaussian(),
                                data = data.frame(mu_pseudo_A_star = mu_pseudo_A_star, A, M, X)[Delta_Y == 1,])
      eta_AXMn_A0 <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, M, X))
      eta_AXMn_A1 <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, M, X))
    }else if(!is.null(SL_eta_AXM)){
      set.seed(seed)
      # since observations with Delta_M=0 will not be considered in the next step
      # so add restrictions here to better depict data with Delta_M=1 # [Delta_M==1]
      eta_AXM_fit <- SuperLearner::SuperLearner(Y = mu_pseudo_A_star[Delta_Y == 1], X = data.frame(A, M, X)[Delta_Y == 1,],
                                                family = gaussian(), 
                                                SL.library = SL_eta_AXM,
                                                method = tmp_method.CC_LS(),
                                                control = list(saveCVFitLibrary = TRUE))
      eta_AXMn_A0 <- stats::predict(eta_AXM_fit, type = "response", newdata = data.frame(A = 0, M, X))[[1]]
      eta_AXMn_A1 <- stats::predict(eta_AXM_fit, type = "response", newdata = data.frame(A = 1, M, X))[[1]]
    }else{
      set.seed(seed)
      eta_AXM_fit <- hal9001::fit_hal(X = data.frame(A, M, X)[Delta_Y == 1,], Y = mu_pseudo_A_star[Delta_Y == 1]) 
      eta_AXMn_A0 <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 0, M, X))
      eta_AXMn_A1 <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 1, M, X))
    }
    
    # fit pseudo-outcome regression for eta_AXZ, conditioning on disease status and baseline covariates
    if(!is.null(glm_formula$xi_AX)){
      xi_fit <- stats::glm(paste0("eta_AXZn_A ~ ", glm_formula$xi_AX), family = gaussian(), 
                          data = data.frame(eta_AXZn_A = eta_AXZn_A, A, X))
      xi_AXn_A0 <- stats::predict(xi_fit, newdata = data.frame(A = 0, X))
      xi_AXn_A1 <- stats::predict(xi_fit, newdata = data.frame(A = 1, X))
    }else{
      set.seed(seed)
      xi_fit <- SuperLearner::SuperLearner(Y = eta_AXZn_A, X = data.frame(A, X),
                                           family = gaussian(),
                                           SL.library = SL_xi_AX,
                                           method = tmp_method.CC_LS(),
                                           control = list(saveCVFitLibrary = TRUE))
      xi_AXn_A0 <- stats::predict(xi_fit, newdata = data.frame(A = 0, X))[[1]]
      xi_AXn_A1 <- stats::predict(xi_fit, newdata = data.frame(A = 1, X))[[1]]
    }
    
    eif_A0 <- make_full_data_eif(a = 0, A = A, Delta_Y = Delta_Y, Delta_M = Delta_M, Y = Y[,j], 
                                 gA = 1 - gAn_1, gDM = gDMn_1_A0, gDY_AXZ = gDYn_1_A0XZ,
                                 eta_AXZ = eta_AXZn_A0, eta_AXM = eta_AXMn_A0, 
                                 xi_AX = xi_AXn_A0, 
                                 mu = mu_MXZn_A0, pMXD = pMXDn_A0, pMXZ = pMXZn_A0)
    
    eif_A1 <- make_full_data_eif(a = 1, A = A, Delta_Y = Delta_Y, Delta_M = Delta_M, Y = Y[,j], 
                                 gA = gAn_1, gDM = gDMn_1_A0, gDY_AXZ = gDYn_1_A1XZ,
                                 eta_AXZ = eta_AXZn_A1, eta_AXM = eta_AXMn_A1, 
                                 xi_AX = xi_AXn_A1, 
                                 mu = mu_MXZn_A1, pMXD = pMXDn_A0, pMXZ = pMXZn_A1)
    
    # one-step estimator
    est[,j] <- c(mean(xi_AXn_A0) + mean(eif_A0), mean(xi_AXn_A1) + mean(eif_A1))
    
    # covariance matrix
    eif_mat[[j]] <- cbind(eif_A0, eif_A1)
    colnames(eif_mat[[j]]) <- c("eif_A0", "eif_A1")
    cov_mat[[j]] <- cov(eif_mat[[j]]) / n
  }
  
  # motion-controlled associations
  adj_association <- est[2, ] - est[1, ]
  
  # output
  out <- list(est = est,
              adj_association = adj_association,
              eif_mat = eif_mat,
              cov_mat = cov_mat
  )
  
  return(out)
}

