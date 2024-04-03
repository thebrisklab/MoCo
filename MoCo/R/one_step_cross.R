#' Non-parametric efficient motion-adjusted functional connectivity estimators
#' 
#' @details Compute non-parametric efficient motion-adjusted functional connectivity estimators based on the proposed method.
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
#' @param HAL_pMX Specifies whether to estimate p(m | a, x) and p(m | a, x, Delta_M=1) using the highly adaptive lasso conditional density estimation method. Defaults to \code{TRUE}. 
#' @param HAL_pMXZ Specifies whether to estimate p(m | a, x, z) and p(m | a, x, z, Delta_M=1) using the highly adaptive lasso conditional density estimation method. Defaults to \code{TRUE}. 
#' 
#' @param HAL_options Additional options for highly adaptive lasso (HAL) method.
#'                   - \code{max_degree}: The highest order of interaction terms for generating basis functions (passed to \code{haldensify}).
#'                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression (passed to \code{haldensify}).
#'                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions (passed to \code{haldensify}).
#' 
#' @param cv_folds a numeric value indicates the number(s) of partitions in cross-fitting. Defaults to 5. 
#' 
#' @importFrom SuperLearner SuperLearner
#' @importFrom haldensify haldensify
#' 
#' @return A list with named entries 
#' \describe{
#'   \item{est}{A two times p matrix showing the one-step estimators of the control group and the disease group for each functional connectivity of interest, respectively.}
#'   \item{adj_association}{A p-length vector showing the motion-adjusted association for each functional connectivity of interest, respectively.}
#'   \item{eif_mat}{A p-length list of the estimated EIF evaluated on the observations of the control and disease group, respectively.}
#'   \item{cov_mat}{A p-length list of the estimated covariance matrix of the one-step estimators.}
#' }
#' 
#' @export

fit_mechanism <- function(
  train_dat, 
  valid_dat,
  SL_library = c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.glm.interaction","SL.ranger", "SL.xgboost","SL.mean"),
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
  HAL_options = list(max_degree = 6, lambda_seq = exp(seq(-1, -10, length = 100)), num_knots = c(1000, 500, 250)),
  seed = 1, 
  ...
){
  # number of edges
  p <- ncol(train_dat$Y)
  
  # motion and motion indicator of training and valid data
  M_train <- Reduce(c, train_dat$M)
  Delta_M_train <- Reduce(c, train_dat$Delta_M)
  M_valid <- Reduce(c, valid_dat$M)
  Delta_M_valid <- Reduce(c, valid_dat$Delta_M)
  Delta_Y_train <- Reduce(c, train_dat$Delta_Y)
  Delta_Y_valid <- Reduce(c, valid_dat$Delta_Y)
  
  # fit regression for propensity score 
  # disease status, conditional on baseline covariates
  if(!is.null(glm_formula$gA)){
    gA_fit <- stats::glm(paste0("A ~ ", glm_formula$gA), family = binomial(), 
                        data = data.frame(train_dat$A, train_dat$X))
    gAn_1 <- stats::predict(gA_fit, type = "response", newdata = valid_dat$X)
  }else{
    set.seed(seed)
    if(ncol(train_dat$X) == 1){
      SL_gA <- SL_library[SL_library != "SL.glmnet"]
    }else{
      SL_gA <- SL_library
    }
    gA_fit <- SuperLearner::SuperLearner(Y = train_dat$A$A, X = train_dat$X,
                                         family = binomial(), 
                                         SL.library = SL_gA,
                                         method = tmp_method.CC_nloglik(),
                                         control = list(saveCVFitLibrary = TRUE))
    gAn_1 <- stats::predict(gA_fit, type = "response", newdata = valid_dat$X)[[1]]
  }
  
  # fit binary mediator regression
  if(!is.null(glm_formula$gDM)){
    gDM_fit <- stats::glm(paste0("Delta_M ~ ", glm_formula$gDM), family = binomial(), data = data.frame(train_dat$Delta_M, train_dat$A, train_dat$X))
    gDMn_1_A0_train <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, train_dat$X))
    gDMn_1_A0_valid <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X))
  }else{
    set.seed(seed)
    gDM_fit <- SuperLearner::SuperLearner(Y = Delta_M_train, X = data.frame(train_dat$A, train_dat$X),
                                          family = binomial(), 
                                          SL.library = SL_library,
                                          method = tmp_method.CC_nloglik(),
                                          control = list(saveCVFitLibrary = TRUE))
    gDMn_1_A0_train <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, train_dat$X))[[1]]
    gDMn_1_A0_valid <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X))[[1]]
  }
  
  # fit missingness indicator regression
  if(sum(Delta_Y_train) == length(Delta_Y_train)){
    gDYn_1_AXZ_train <- gDYn_1_A0XZ_train <- gDYn_1_A1XZ_train <- rep(1, length(Delta_Y_train))
    gDYn_1_AXZ_valid <- gDYn_1_A0XZ_valid <- gDYn_1_A1XZ_valid <- rep(1, length(Delta_Y_valid))
    
    gDYn_1_AX_train <- rep(1, length(Delta_Y_train))
    gDYn_1_AX_valid <- rep(1, length(Delta_Y_valid))
  }else{
    # probability of non-missing conditioning on disease status A and baseline covariates X 
    if(!is.null(glm_formula$gDY_AX)){
      gDY_AX_fit <- stats::glm(paste0("Delta_Y ~ ", glm_formula$gDY_AX), family = binomial(), 
                               data = data.frame(train_dat$Delta_Y, train_dat$A, train_dat$X))
      gDYn_1_AX_train <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X))
      gDYn_1_AX_valid <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X))
    }else{
      set.seed(seed)
      gDY_AX_fit <- SuperLearner::SuperLearner(Y = Delta_Y_train, X = data.frame(train_dat$A, train_dat$X),
                                               family = binomial(), 
                                               SL.library = SL_library,
                                               method = tmp_method.CC_nloglik(),
                                               control = list(saveCVFitLibrary = TRUE))
      gDYn_1_AX_train <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X))[[1]]
      gDYn_1_AX_valid <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X))[[1]]
    } 
    
    # probability of non-missing conditioning on disease status A, baseline covariates X and diagnosis-related covariates Z
    if(!is.null(glm_formula$gDY_AXZ)){
      gDY_AXZ_fit <- stats::glm(paste0("Delta_Y ~ ", glm_formula$gDY_AXZ), family = binomial(), 
                                data = data.frame(train_dat$Delta_Y, train_dat$A, train_dat$X, train_dat$Z))
      gDYn_1_AXZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z))
      gDYn_1_A0XZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, train_dat$X, train_dat$Z))
      gDYn_1_A1XZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, train_dat$X, train_dat$Z))
      gDYn_1_AXZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z))
      gDYn_1_A0XZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z))
      gDYn_1_A1XZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z))
    }else{
      set.seed(seed)
      gDY_AXZ_fit <- SuperLearner::SuperLearner(Y = Delta_Y_train, X = data.frame(train_dat$A, train_dat$X, train_dat$Z),
                                                family = binomial(), 
                                                SL.library = SL_library,
                                                method = tmp_method.CC_nloglik(),
                                                control = list(saveCVFitLibrary = TRUE))
      gDYn_1_AXZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z))[[1]]
      gDYn_1_A0XZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, train_dat$X, train_dat$Z))[[1]]
      gDYn_1_A1XZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, train_dat$X, train_dat$Z))[[1]]
      gDYn_1_AXZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z))[[1]]
      gDYn_1_A0XZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z))[[1]]
      gDYn_1_A1XZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z))[[1]]
    }
  }
  
  if(HAL_pMX){
    # density estimation for mediator, conditioning on diease status and baseline covariates p(m|a,x)
    # estimate p(m|a,x) use highly adaptive lasso conditional density estimation method
    # use default n_bins range, use cv to choose number of bins
    pMX_fit <- haldensify::haldensify(
      A = M_train[Delta_Y_train == 1],
      W = data.frame(train_dat$A, train_dat$X)[Delta_Y_train == 1,],  
      max_degree = HAL_options$max_degree,
      lambda_seq = HAL_options$lambda_seq,
      num_knots = HAL_options$num_knots
    ) 
    pMXn_A_train <- rep(NA, length(Delta_M_train))
    pMXn_A_valid <- rep(NA, length(Delta_M_valid))
    pMXn_A_train[Delta_Y_train == 1] <- stats::predict(pMX_fit, new_A = M_train[Delta_Y_train == 1], new_W = data.frame(train_dat$A, train_dat$X)[Delta_Y_train == 1,])
    pMXn_A_valid[Delta_Y_valid == 1] <- stats::predict(pMX_fit, new_A = M_valid[Delta_Y_valid == 1], new_W = data.frame(valid_dat$A, valid_dat$X)[Delta_Y_valid == 1,])
    
    # density estimation for mediator, conditioning on diease status and baseline covariates p(m|0,x,Delta=1)
    # use default n_bins range, use cv to choose number of bins
    pMXD_fit <- haldensify::haldensify(
      A = M_train[Delta_M_train == 1],
      W = data.frame(train_dat$A, train_dat$X)[Delta_M_train == 1,],  
      max_degree = HAL_options$max_degree,
      lambda_seq = HAL_options$lambda_seq,
      num_knots = HAL_options$num_knots
    ) 
    pMXDn_A0_train <- rep(NA, length(Delta_M_train))
    pMXDn_A0_train[which(!is.na(M_train))] <- stats::predict(pMXD_fit, new_A = M_train[which(!is.na(M_train))], new_W = data.frame(A = 0, train_dat$X)[which(!is.na(M_train)),], trim_min = 0)
    pMXDn_A0_valid <- rep(NA, length(Delta_M_valid))
    pMXDn_A0_valid[which(!is.na(M_valid))] <- stats::predict(pMXD_fit, new_A = M_valid[which(!is.na(M_valid))], new_W = data.frame(A = 0, valid_dat$X)[which(!is.na(M_valid)),], trim_min = 0)
  }else{
    pMX_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X)[Delta_Y_train == 1,])
    pMXn_A_train <- rep(NA, length(Delta_M_train))
    pMXn_A_valid <- rep(NA, length(Delta_M_valid))
    pMXn_A_train[Delta_Y_train == 1] <- (1/M_train)[Delta_Y_train == 1] * dnorm(log(M_train)[Delta_Y_train == 1], mean = stats::predict(pMX_fit, newdata = data.frame(train_dat$A, train_dat$X)[Delta_Y_train == 1,]), sd = sd(pMX_fit$residuals)) 
    pMXn_A_valid[Delta_Y_valid == 1] <- (1/M_valid)[Delta_Y_valid == 1] * dnorm(log(M_valid)[Delta_Y_valid == 1], mean = stats::predict(pMX_fit, newdata = data.frame(valid_dat$A, valid_dat$X)[Delta_Y_valid == 1,]), sd = sd(pMX_fit$residuals))
    
    pMXD_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X)[Delta_M_train == 1,])
    pMXDn_A0_train <- rep(NA, length(Delta_M_train))
    pMXDn_A0_train[Delta_M_train==0] <- 0
    pMXDn_A0_train[Delta_M_train==1] <- (1/M_train)[Delta_M_train == 1] * dnorm(log(M_train)[Delta_M_train == 1], mean = stats::predict(pMXD_fit, newdata = data.frame(A=0, train_dat$X)[Delta_M_train == 1,]), sd = sd(pMXD_fit$residuals))
    pMXDn_A0_valid <- rep(NA, length(Delta_M_valid))
    pMXDn_A0_valid[Delta_M_valid==0] <- 0
    pMXDn_A0_valid[Delta_M_valid==1] <- (1/M_valid)[Delta_M_valid == 1] * dnorm(log(M_valid)[Delta_M_valid == 1], mean = stats::predict(pMXD_fit, newdata = data.frame(A=0, valid_dat$X)[Delta_M_valid == 1,]), sd = sd(pMXD_fit$residuals))
  }
  
  if(HAL_pMXZ){
    # density estimation for mediator, conditioning on disease status, binary mediator, 
    # baseline covariates, and mediator-outcome confounder p(m|0,x,z)
    # use default n_bins range, use cv to choose number of bins
    pMXZ_fit <- haldensify::haldensify(
      A = M_train[Delta_Y_train == 1],
      W = data.frame(train_dat$A, train_dat$X, train_dat$Z)[Delta_Y_train == 1,], 
      max_degree = HAL_options$max_degree,
      lambda_seq = HAL_options$lambda_seq,
      num_knots = HAL_options$num_knots
    )
    pMXZn_A0_train <- pMXZn_A1_train <- pMXZn_A_train <- rep(NA, length(Delta_M_train))
    pMXZn_A_train[Delta_Y_train==1] <- stats::predict(pMXZ_fit, new_A = M_train[Delta_Y_train==1], new_W = data.frame(train_dat$A, train_dat$X, train_dat$Z)[Delta_Y_train==1,])
    pMXZn_A0_train[Delta_Y_train==1] <- stats::predict(pMXZ_fit, new_A = M_train[Delta_Y_train==1], new_W = data.frame(A = 0, train_dat$X, train_dat$Z)[Delta_Y_train==1,])
    pMXZn_A1_train[Delta_Y_train==1] <- stats::predict(pMXZ_fit, new_A = M_train[Delta_Y_train==1], new_W = data.frame(A = 1, train_dat$X, train_dat$Z)[Delta_Y_train==1,])
    pMXZn_A0_valid <- pMXZn_A1_valid <- pMXZn_A_valid <- rep(NA, length(Delta_M_valid))
    pMXZn_A_valid[Delta_Y_valid==1] <- stats::predict(pMXZ_fit, new_A = M_valid[Delta_Y_valid==1], new_W = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z)[Delta_Y_valid==1,])
    pMXZn_A0_valid[Delta_Y_valid==1] <- stats::predict(pMXZ_fit, new_A = M_valid[Delta_Y_valid==1], new_W = data.frame(A = 0, valid_dat$X, valid_dat$Z)[Delta_Y_valid==1,])
    pMXZn_A1_valid[Delta_Y_valid==1] <- stats::predict(pMXZ_fit, new_A = M_valid[Delta_Y_valid==1], new_W = data.frame(A = 1, valid_dat$X, valid_dat$Z)[Delta_Y_valid==1,])
    
    # density estimation for mediator, conditioning on disease status, binary mediator, 
    # baseline covariates, and mediator-outcome confounder p(m|0,x,z,Delta=1)
    pMXZD_fit <- haldensify::haldensify(
      A = M_train[Delta_M_train == 1],
      W = data.frame(train_dat$A, train_dat$X, train_dat$Z)[Delta_M_train == 1,], 
      max_degree = HAL_options$max_degree,
      lambda_seq = HAL_options$lambda_seq,
      num_knots = HAL_options$num_knots
    )
    pMXZDn_A_train <- rep(NA, length(Delta_M_train))
    pMXZDn_A_train[which(!is.na(M_train))] <- stats::predict(pMXZD_fit, new_A = M_train[which(!is.na(M_train))], new_W = data.frame(train_dat$A, train_dat$X, train_dat$Z)[which(!is.na(M_train)),], trim_min = 0)
    pMXZDn_A_valid <- rep(NA, length(Delta_M_valid))
    pMXZDn_A_valid[which(!is.na(M_valid))] <- stats::predict(pMXZD_fit, new_A = M_valid[which(!is.na(M_valid))], new_W = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z)[which(!is.na(M_valid)),], trim_min = 0)
  }else{
    pMXZ_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMXZ), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X, train_dat$Z)[Delta_Y_train == 1,])
    pMXZn_A0_train <- pMXZn_A1_train <- pMXZn_A_train <- rep(NA, length(Delta_M_train))
    pMXZn_A_train[Delta_Y_train==1] <- (1/M_train)[Delta_Y_train==1] * dnorm(log(M_train)[Delta_Y_train==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z)[Delta_Y_train==1, ]), sd = sd(pMXZ_fit$residuals))
    pMXZn_A0_train[Delta_Y_train==1] <- (1/M_train)[Delta_Y_train==1] * dnorm(log(M_train)[Delta_Y_train==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 0, train_dat$X, train_dat$Z)[Delta_Y_train==1, ]), sd = sd(pMXZ_fit$residuals))
    pMXZn_A1_train[Delta_Y_train==1] <- (1/M_train)[Delta_Y_train==1] * dnorm(log(M_train)[Delta_Y_train==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 1, train_dat$X, train_dat$Z)[Delta_Y_train==1, ]), sd = sd(pMXZ_fit$residuals))      
    pMXZn_A0_valid <- pMXZn_A1_valid <- pMXZn_A_valid <- rep(NA, length(Delta_M_valid))
    pMXZn_A_valid[Delta_Y_valid==1] <- (1/M_valid)[Delta_Y_valid==1] * dnorm(log(M_valid)[Delta_Y_valid==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z)[Delta_Y_valid==1, ]), sd = sd(pMXZ_fit$residuals))
    pMXZn_A0_valid[Delta_Y_valid==1] <- (1/M_valid)[Delta_Y_valid==1] * dnorm(log(M_valid)[Delta_Y_valid==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z)[Delta_Y_valid==1, ]), sd = sd(pMXZ_fit$residuals))
    pMXZn_A1_valid[Delta_Y_valid==1] <- (1/M_valid)[Delta_Y_valid==1] * dnorm(log(M_valid)[Delta_Y_valid==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z)[Delta_Y_valid==1, ]), sd = sd(pMXZ_fit$residuals))      
    
    pMXZD_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMXZ), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X, train_dat$Z)[Delta_M_train == 1,])
    pMXZDn_A_train <- rep(NA, length(Delta_M_train))
    pMXZDn_A_train[Delta_M_train==0] <- 0
    pMXZDn_A_train[Delta_M_train==1] <- (1/M_train)[Delta_M_train==1] * dnorm(log(M_train)[Delta_M_train==1], mean = stats::predict(pMXZD_fit, newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z)[Delta_M_train==1,]), 
                                              sd = sd(pMXZD_fit$residuals))
    pMXZDn_A_valid <- rep(NA, length(Delta_M_valid))
    pMXZDn_A_valid[Delta_M_valid==0] <- 0
    pMXZDn_A_valid[Delta_M_valid==1] <- (1/M_valid)[Delta_M_valid==1] * dnorm(log(M_valid)[Delta_M_valid==1], mean = stats::predict(pMXZD_fit, newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z)[Delta_M_valid==1,]), 
                                              sd = sd(pMXZD_fit$residuals))
  }
  
  # define parameters for storing results
  est <- matrix(nrow = 2, ncol = p)
  eif_mat <- list(p)
  cov_mat <- list(p)
  for(j in 1:p){
    # fit outcome regression
    if(!is.null(glm_formula$mu_AMXZ)){
      mu_AMXZ_fit <- stats::glm(paste0("Y ~ ", glm_formula$mu_AMXZ), family = gaussian(), data = data.frame(Y = train_dat$Y[,j], train_dat$A, train_dat$M, train_dat$X, train_dat$Z)[Delta_Y_train==1,])
      mu_MXZn_A0_train <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, train_dat$M, train_dat$X, train_dat$Z))
      mu_MXZn_A1_train <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, train_dat$M, train_dat$X, train_dat$Z))
      mu_MXZn_A0_valid <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, valid_dat$M, valid_dat$X, valid_dat$Z))
      mu_MXZn_A1_valid <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, valid_dat$M, valid_dat$X, valid_dat$Z))
      mu_MXZn_A_train <- stats::predict(mu_AMXZ_fit, newdata = data.frame(train_dat$A, train_dat$M, train_dat$X, train_dat$Z))
      mu_MXZn_A_valid <- stats::predict(mu_AMXZ_fit, newdata = data.frame(valid_dat$A, valid_dat$M, valid_dat$X, valid_dat$Z))
    }else{
      set.seed(seed)
      mu_AMXZ_fit <- SuperLearner::SuperLearner(Y = Reduce(c, train_dat$Y[,j])[Delta_Y_train==1], X = data.frame(train_dat$A, train_dat$M, train_dat$X, train_dat$Z)[Delta_Y_train==1,],
                                                family = gaussian(), 
                                                SL.library = SL_library,
                                                method = tmp_method.CC_LS(),
                                                control = list(saveCVFitLibrary = TRUE))
      mu_MXZn_A0_train <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, train_dat$M, train_dat$X, train_dat$Z))[[1]]
      mu_MXZn_A1_train <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, train_dat$M, train_dat$X, train_dat$Z))[[1]]
      mu_MXZn_A0_valid <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, valid_dat$M, valid_dat$X, valid_dat$Z))[[1]]
      mu_MXZn_A1_valid <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, valid_dat$M, valid_dat$X, valid_dat$Z))[[1]]
      mu_MXZn_A_train <- stats::predict(mu_AMXZ_fit, newdata = data.frame(train_dat$A, train_dat$M, train_dat$X, train_dat$Z))[[1]]
      mu_MXZn_A_valid <- stats::predict(mu_AMXZ_fit, newdata = data.frame(valid_dat$A, valid_dat$M, valid_dat$X, valid_dat$Z))[[1]]
    }
    
    # fit additional pseudo-outcome regression for QY*qM/qMZ, conditioning on disease status, binary mediator and baseline covariates
    mu_pseudo_A_train <- mu_MXZn_A_train*pMXDn_A0_train/pMXZDn_A_train
    if(!is.null(glm_formula$eta_AXZ)){
      eta_AXZ_fit <- stats::glm(paste0("mu_pseudo_A ~ ", glm_formula$eta_AXZ), family = gaussian(),
                                data = data.frame(mu_pseudo_A = mu_pseudo_A_train, train_dat$A, train_dat$X, train_dat$Z)[Delta_M_train==1,])
      eta_AXZn_A_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z))
      eta_AXZn_A0_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, train_dat$X, train_dat$Z))
      eta_AXZn_A1_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, train_dat$X, train_dat$Z))
      
      eta_AXZn_A_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z))
      eta_AXZn_A0_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z))
      eta_AXZn_A1_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z))
    }else{
      set.seed(seed)
      eta_AXZ_fit <- SuperLearner::SuperLearner(Y = mu_pseudo_A_train[Delta_M_train==1],
                                                X = data.frame(train_dat$A, train_dat$X, train_dat$Z)[Delta_M_train==1,],
                                                family = gaussian(), 
                                                SL.library = SL_library,
                                                method = tmp_method.CC_LS(),
                                                control = list(saveCVFitLibrary = TRUE))
      eta_AXZn_A_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z))[[1]]
      eta_AXZn_A0_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, train_dat$X, train_dat$Z))[[1]]
      eta_AXZn_A1_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, train_dat$X, train_dat$Z))[[1]]
      
      eta_AXZn_A_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z))[[1]]
      eta_AXZn_A0_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z))[[1]]
      eta_AXZn_A1_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z))[[1]]
    }
    
    # fit pseudo-outcome regression for Qd, conditioning on disease status, mediator and baseline covariates
    mu_pseudo_A_star_train <- mu_MXZn_A_train*pMXn_A_train/pMXZn_A_train*gDYn_1_AX_train/gDYn_1_AXZ_train
    # fit the regression
    if(!is.null(glm_formula$eta_AXM)){
      eta_AXM_fit <- stats::glm(paste0("mu_pseudo_A_star ~ ", glm_formula$eta_AXM), family = gaussian(),
                                data = data.frame(mu_pseudo_A_star = mu_pseudo_A_star_train, train_dat$A, train_dat$M, train_dat$X)[Delta_Y_train==1,])
      eta_AXMn_A0_train <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, train_dat$M, train_dat$X))
      eta_AXMn_A1_train <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, train_dat$M, train_dat$X))
      
      eta_AXMn_A0_valid <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, valid_dat$M, valid_dat$X))
      eta_AXMn_A1_valid <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, valid_dat$M, valid_dat$X))
    }else if(!is.null(SL_library)){
      set.seed(seed)
      # since observations with Delta=0 will not be considered in the next step
      # so add restrictions here to better depict data with Delta=1
      eta_AXM_fit <- SuperLearner::SuperLearner(Y = mu_pseudo_A_star_train[Delta_Y_train==1], 
                                                X = data.frame(train_dat$A, train_dat$M, train_dat$X)[Delta_Y_train==1,],
                                                family = gaussian(),
                                                SL.library = SL_library,
                                                method = tmp_method.CC_LS(),
                                                control = list(saveCVFitLibrary = TRUE))
      eta_AXMn_A0_train <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, train_dat$M, train_dat$X))[[1]]
      eta_AXMn_A1_train <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, train_dat$M, train_dat$X))[[1]]
      
      eta_AXMn_A0_valid <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, valid_dat$M, valid_dat$X))[[1]]
      eta_AXMn_A1_valid <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, valid_dat$M, valid_dat$X))[[1]]
    }else{
      set.seed(seed)
      eta_AXM_fit <- hal9001::fit_hal(X = data.frame(train_dat$A, train_dat$M, train_dat$X)[Delta_Y_train==1,], Y = mu_pseudo_A0_star[Delta_Y_train==1]) 
      eta_AXMn_A0_train <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 0, train_dat$M, train_dat$X))
      eta_AXMn_A1_train <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 1, train_dat$M, train_dat$X))
      
      eta_AXMn_A0_valid <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 0, valid_dat$M, valid_dat$X))
      eta_AXMn_A1_valid <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 1, valid_dat$M, valid_dat$X))
    }
    
    # fit pseudo-outcome regression for u_star, conditioning on disease status and baseline covariates
    if(!is.null(glm_formula$xi_AX)){
      xi_fit <- stats::glm(paste0("eta_AXZn_A ~ ", glm_formula$xi_AX), family = gaussian(),
                           data = data.frame(eta_AXZn_A = eta_AXZn_A_train, train_dat$A, train_dat$X))
      xi_AXn_A0 <- stats::predict(xi_fit, newdata = data.frame(A = 0, valid_dat$X))
      xi_AXn_A1 <- stats::predict(xi_fit, newdata = data.frame(A = 1, valid_dat$X))
    }else{
      set.seed(seed)
      xi_fit <- SuperLearner::SuperLearner(Y = eta_AXZn_A_train, X = data.frame(train_dat$A, train_dat$X),
                                           family = gaussian(), 
                                           SL.library = SL_library,
                                           method = tmp_method.CC_LS(),
                                           control = list(saveCVFitLibrary = TRUE))
      xi_AXn_A0 <- stats::predict(xi_fit, newdata = data.frame(A = 0, valid_dat$X))[[1]]
      xi_AXn_A1 <- stats::predict(xi_fit, newdata = data.frame(A = 1, valid_dat$X))[[1]]
    }
    
    # calculate efficient influence function
    eif_A0 <- make_full_data_eif(a = 0, A = valid_dat$A$A, Delta_Y = Delta_Y_valid, Delta_M = Delta_M_valid, Y = valid_dat$Y[,j], 
                                 gA = 1 - gAn_1, gDM = gDMn_1_A0_valid, gDY_AXZ = gDYn_1_A0XZ_valid,
                                 eta_AXZ = eta_AXZn_A0_valid, eta_AXM = eta_AXMn_A0_valid, 
                                 xi_AX = xi_AXn_A0, 
                                 mu = mu_MXZn_A0_valid, pMXD = pMXDn_A0_valid, pMXZ = pMXZn_A0_valid)
    
    eif_A1 <- make_full_data_eif(a = 1, A = valid_dat$A$A, Delta_Y = Delta_Y_valid, Delta_M = Delta_M_valid, Y = valid_dat$Y[,j], 
                                 gA = gAn_1, gDM = gDMn_1_A0_valid, gDY_AXZ = gDYn_1_A1XZ_valid,
                                 eta_AXZ = eta_AXZn_A1_valid, eta_AXM = eta_AXMn_A1_valid, 
                                 xi_AX = xi_AXn_A1, 
                                 mu = mu_MXZn_A1_valid, pMXD = pMXDn_A0_valid, pMXZ = pMXZn_A1_valid)
    
    # one-step estimator
    est[,j] <- c(mean(xi_AXn_A0) + mean(eif_A0), mean(xi_AXn_A1) + mean(eif_A1))
    
    # covariance matrix
    eif_mat[[j]] <- cbind(eif_A0, eif_A1)
    colnames(eif_mat[[j]]) <- c("eif_A0", "eif_A1")
  }


  # output
  out <- list(est = est,
             eif_mat = eif_mat
  )
  
  return(out)
}

#' Non-parametric efficient motion-adjusted functional connectivity estimators
#' 
#' @details Compute non-parametric efficient motion-adjusted functional connectivity estimators based on the proposed method.
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
#' @param SL_library SuperLearner library for estimating nuisance regressions. Defaults to SL_library if not specified.
#' @param glm_formula All glm formulas default to NULL, indicating SuperLearner will be used for nuisance regressions.
#'                    - \code{gA}: GLM formula for estimating the propensity score.
#'                    - \code{gD}: GLM formula for estimating the probability P(Delta_M = 1 | A, X).
#'                    - \code{mu_AMXZ}: GLM formula for estimating the outcome regression E(Y | A, M, X, Z).
#'                    - \code{xi_AX}: GLM formula for estimating E(eta_AXZ | A=a, X).
#'                    - \code{eta_AXM}: GLM formula for estimating E(Q * p / r | A=a, M, X).
#'                    - \code{pMX}: GLM formula for estimating p(m | a, x) and p(m | a, x, Delta_M = 1), assuming M follows a normal distribution.
#'                    - \code{pMXZ}: GLM formula for estimating p(m | a, x, z) and p(m | a, x, z, Delta_M = 1), assuming M follows a normal distribution.
#'                    
#' @param HAL_pMX Specifies whether to estimate p(m | a, x) and p(m | a, x, Delta_M=1) using the highly adaptive lasso conditional density estimation method. Defaults to \code{TRUE}. 
#' @param HAL_pMXZ Specifies whether to estimate p(m | a, x, z) and p(m | a, x, z, Delta_M=1) using the highly adaptive lasso conditional density estimation method. Defaults to \code{TRUE}. 
#' 
#' @param HAL_options Additional options for highly adaptive lasso (HAL) method.
#'                   - \code{max_degree}: The highest order of interaction terms for generating basis functions (passed to \code{haldensify}).
#'                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression (passed to \code{haldensify}).
#'                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions (passed to \code{haldensify}).
#' 
#' @param cv_folds a numeric value indicates the number(s) of partitions in cross-fitting. Defaults to 5. 
#' 
#' @importFrom SuperLearner SuperLearner
#' @importFrom haldensify haldensify
#' 
#' @return A list with named entries 
#' \describe{
#'   \item{est}{A two times p matrix showing the one-step estimators of the control group and the disease group for each functional connectivity of interest, respectively.}
#'   \item{adj_association}{A p-length vector showing the motion-adjusted association for each functional connectivity of interest, respectively.}
#'   \item{eif_mat}{A p-length list of the estimated EIF evaluated on the observations of the control and disease group, respectively.}
#'   \item{cov_mat}{A p-length list of the estimated covariance matrix of the one-step estimators.}
#' }

one_step_cross <- function(
    X, Z, A, M, Y, 
    Delta_M, 
    thresh = NULL,
    Delta_Y,
    SL_library = c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.glm.interaction","SL.ranger", "SL.xgboost","SL.mean"),
    glm_formula = list(gA = NULL, 
                       gDM = NULL,
                       gDY_X = NULL,
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
    seed = 1, 
    cv_folds = 5,
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
  
  # the original dataset
  dat <- list(X = data.frame(X), Z = data.frame(Z), A = data.frame(A), 
              M = data.frame(M), Delta_M = data.frame(Delta_M), 
              Delta_Y = data.frame(Delta_Y), Y = data.frame(Y))
  
  # divide the dataset into train_dat and valid_dat
  set.seed(seed)
  folds <- caret::createFolds(1:n, k = cv_folds, list = TRUE, returnTrain = FALSE)
  
  # create empty matrix to store results
  est <- matrix(0, nrow = 2, ncol = p)
  # create empty vector to store eifs
  eif_mat <- list(p)
  for(j in 1:p){
    eif_mat[[j]] <- matrix(nrow = n, ncol = 2)
  }
  
  # cross fitting
  for(i in 1:cv_folds){
    rslt_cross <- fit_mechanism(
      train_dat = lapply(dat, function(x){x[-folds[[i]], , drop = FALSE]}), 
      valid_dat = lapply(dat, function(x){x[folds[[i]], , drop = FALSE]}),
      SL_library = SL_library,
      glm_formula = glm_formula,
      HAL_pMX = HAL_pMX,
      HAL_pMXZ = HAL_pMXZ,
      HAL_options = HAL_options,
      seed = seed
    )

    # store results
    est <- est + rslt_cross$est
    for(j in 1:p){
      eif_mat[[j]][folds[[i]],] <- rslt_cross$eif_mat[[j]]
    }
  }
  
  # calculate one-step estimator
  est <- est/cv_folds 
  
  # motion-adjusted associations
  adj_association <- est[2, ] - est[1, ]
  
  # calculate covariance matrix
  cov_mat <- list(p)
  for(j in 1:p){
    colnames(eif_mat[[j]]) <- c('est_A0', 'est_A1')
    cov_mat[[j]] <- cov(eif_mat[[j]])/n
  }
  
  # output
  out <- list(est = est,
              adj_association = adj_association,
              eif_mat = eif_mat,
              cov_mat = cov_mat
  )
}

