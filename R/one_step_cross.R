#' Non-parametric efficient motion-controlled brain phenotype estimators
#' 
#' @details Compute non-parametric efficient motion-controlled brain phenotype estimators (MoCo) using one-step cross-fitting
#' 
#' @param train_dat A \code{list} containing the observed data used for training, with items in the order (A, X, Z, M, Y, Delta). 
#' @param valid_dat A \code{list} containing the validation data used for prediction, with items in the order (A, X, Z, M, Y, Delta). 
#' @param M_indicator A binary value indicating whether motion value is available.
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
#'                    - \code{pMX}: GLM formula for estimating p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M = 1), assuming M follows a normal distribution.
#'                    - \code{pMXZ}: GLM formula for estimating p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M = 1), assuming M follows a normal distribution.
#'                    
#' @param HAL_pMX Specifies whether to estimate p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M=1) using the highly adaptive lasso conditional density estimation method. 
#' Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMX} option in \code{glm_formula}, such as \code{pMX = "."}.
#' @param HAL_pMXZ Specifies whether to estimate p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M=1) using the highly adaptive lasso conditional density estimation method. 
#' Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMXZ} option in \code{glm_formula}, such as \code{pMXZ = "."}.
#'  
#' @param HAL_options Additional options for highly adaptive lasso (HAL) method.
#'                   - \code{max_degree}: The highest order of interaction terms for generating basis functions (passed to \code{haldensify}).
#'                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression (passed to \code{haldensify}).
#'                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions (passed to \code{haldensify}).
#' 
#' 
#' @param density_method,Delta_M_threshold,gamlss_options Motion density method, the threshold that defines \code{Delta_M}, and GAMLSS options. See \code{\link{moco}}.
#' @param seed Random seed set before each SuperLearner fit (one value of \code{seed_rng} in \code{moco()}).
#' @param ... Not used.
#' 
#' @import SuperLearner 
#' @import haldensify 
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
    M_indicator,
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
      max_degree = 3, 
      lambda_seq = exp(seq(-1, -10, length = 100)), 
      num_knots = c(1000, 500, 250)
    ),
    density_method = NULL,
    Delta_M_threshold = NULL,
    gamlss_options = NULL,
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
  
  # number of samples
  n_train <- nrow(train_dat$Y)
  n_valid <- nrow(valid_dat$Y)
  
  # non-missing motion index
  index_nomissingM_train <- which(!is.na(M_train))
  index_nomissingM_valid <- which(!is.na(M_valid))
  
  # motion density methods for pMX/pMXD and pMXZ/pMXZD: "hal", "glm" or "gamlss"
  dens_method <- resolve_density_method(density_method, HAL_pMX, HAL_pMXZ)
  if(dens_method["pMX"] == "glm" && is.null(glm_formula$pMX)){glm_formula$pMX <- "."}
  if(dens_method["pMXZ"] == "glm" && is.null(glm_formula$pMXZ)){glm_formula$pMXZ <- "."}
  if(any(dens_method == "gamlss")){
    gg_check_packages()
    gg_opts <- gg_check_options(gamlss_options)
  }
  density_log <- NULL
  weight_diagnostics <- NULL
  
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
                        data = data.frame(train_dat$A, train_dat$X))
    gAn_1 <- stats::predict(gA_fit, type = "response", newdata = valid_dat$X)
  }else{
    set.seed(seed)
    if(ncol(train_dat$X) == 1){
      SL_gA <- SL_gA[SL_gA != "SL.glmnet"]
    }
    gA_fit <- SuperLearner::SuperLearner(Y = train_dat$A$A, X = train_dat$X,
                                         family = binomial(), 
                                         SL.library = SL_gA,
                                         method = tmp_method.CC_nloglik(),
                                         control = list(saveCVFitLibrary = TRUE))
    gAn_1 <- stats::predict(gA_fit, type = "response", newdata = valid_dat$X, onlySL = TRUE)[[1]]
  }
  
  # fit missingness indicator regression
  if(sum(Delta_Y_train) == n_train){
    gDYn_1_AXZ_train <- gDYn_1_A0XZ_train <- gDYn_1_A1XZ_train <- rep(1, n_train)
    gDYn_1_AXZ_valid <- gDYn_1_A0XZ_valid <- gDYn_1_A1XZ_valid <- rep(1, n_valid)
  }else{
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
                                                SL.library = SL_gDY_AXZ,
                                                method = tmp_method.CC_nloglik(),
                                                control = list(saveCVFitLibrary = TRUE))
      gDYn_1_AXZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z), onlySL = TRUE)[[1]]
      gDYn_1_A0XZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, train_dat$X, train_dat$Z), onlySL = TRUE)[[1]]
      gDYn_1_A1XZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, train_dat$X, train_dat$Z), onlySL = TRUE)[[1]]
      gDYn_1_AXZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z), onlySL = TRUE)[[1]]
      gDYn_1_A0XZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z), onlySL = TRUE)[[1]]
      gDYn_1_A1XZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z), onlySL = TRUE)[[1]]
    }
  }
  
  if(M_indicator){
    # fit binary mediator regression
    if(!is.null(glm_formula$gDM)){
      gDM_fit <- stats::glm(paste0("Delta_M ~ ", glm_formula$gDM), family = binomial(), data = data.frame(train_dat$Delta_M, train_dat$A, train_dat$X))
      gDMn_1_A0_train <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, train_dat$X))
      gDMn_1_A0_valid <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X))
    }else{
      set.seed(seed)
      gDM_fit <- SuperLearner::SuperLearner(Y = Delta_M_train, X = data.frame(train_dat$A, train_dat$X),
                                            family = binomial(), 
                                            SL.library = SL_gDM,
                                            method = tmp_method.CC_nloglik(),
                                            control = list(saveCVFitLibrary = TRUE))
      gDMn_1_A0_train <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, train_dat$X), onlySL = TRUE)[[1]]
      gDMn_1_A0_valid <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X), onlySL = TRUE)[[1]]
    }
    
    # fit missingness indicator regression
    if(sum(Delta_Y_train) == n_train){
      gDYn_1_AX_train <- rep(1, n_train)
      gDYn_1_AX_valid <- rep(1, n_valid)
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
                                                 SL.library = SL_gDY_AX,
                                                 method = tmp_method.CC_nloglik(),
                                                 control = list(saveCVFitLibrary = TRUE))
        gDYn_1_AX_train <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X), onlySL = TRUE)[[1]]
        gDYn_1_AX_valid <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X), onlySL = TRUE)[[1]]
      } 
    }
    
    if(dens_method["pMX"] == "hal"){
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
      pMXn_A_train <- rep(NA, n_train)
      pMXn_A_valid <- rep(NA, n_valid)
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
      pMXDn_A0_train <- rep(NA, n_train)
      pMXDn_A0_train[index_nomissingM_train] <- stats::predict(pMXD_fit, new_A = M_train[index_nomissingM_train], new_W = data.frame(A = 0, train_dat$X)[index_nomissingM_train,], trim_min = 0)
      pMXDn_A0_valid <- rep(NA, n_valid)
      pMXDn_A0_valid[index_nomissingM_valid] <- stats::predict(pMXD_fit, new_A = M_valid[index_nomissingM_valid], new_W = data.frame(A = 0, valid_dat$X)[index_nomissingM_valid,], trim_min = 0)
    }else if(dens_method["pMX"] == "gamlss"){
      # GAMLSS generalized gamma densities p(m|a,x,Delta_Y=1) and p(m|0,x,Delta_M=1);
      # pMXD is right-truncated at Delta_M_threshold when it is given
      W_train <- data.frame(train_dat$A, train_dat$X)
      W_valid <- data.frame(valid_dat$A, valid_dat$X)
      it <- index_nomissingM_train
      iv <- index_nomissingM_valid
      pMX_fit <- gg_component(M_train[Delta_Y_train == 1], W_train[Delta_Y_train == 1, , drop = FALSE], "pMX", gg_opts)
      pMXn_A_train <- rep(NA, n_train)
      pMXn_A_valid <- rep(NA, n_valid)
      pMXn_A_train[Delta_Y_train == 1] <- gg_density(pMX_fit, M_train[Delta_Y_train == 1], W_train[Delta_Y_train == 1, , drop = FALSE], gg_opts)
      pMXn_A_valid[Delta_Y_valid == 1] <- gg_density(pMX_fit, M_valid[Delta_Y_valid == 1], W_valid[Delta_Y_valid == 1, , drop = FALSE], gg_opts)
      
      pMXD_fit <- gg_component(M_train[Delta_M_train == 1], W_train[Delta_M_train == 1, , drop = FALSE], "pMXD", gg_opts, threshold = Delta_M_threshold)
      pMXDn_A0_train <- rep(NA, n_train)
      pMXDn_A0_train[it] <- gg_density(pMXD_fit, M_train[it], data.frame(A = 0, train_dat$X)[it, , drop = FALSE], gg_opts, zero_outside = TRUE)
      pMXDn_A0_valid <- rep(NA, n_valid)
      pMXDn_A0_valid[iv] <- gg_density(pMXD_fit, M_valid[iv], data.frame(A = 0, valid_dat$X)[iv, , drop = FALSE], gg_opts, zero_outside = TRUE)
      density_log <- rbind(density_log, pMX_fit$log, pMXD_fit$log)
    }else{
      pMX_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X)[Delta_Y_train == 1,])
      pMXn_A_train <- rep(NA, n_train)
      pMXn_A_valid <- rep(NA, n_valid)
      pMXn_A_train[Delta_Y_train == 1] <- (1/M_train)[Delta_Y_train == 1] * dnorm(log(M_train)[Delta_Y_train == 1], mean = stats::predict(pMX_fit, newdata = data.frame(train_dat$A, train_dat$X)[Delta_Y_train == 1,]), sd = sd(pMX_fit$residuals)) 
      pMXn_A_valid[Delta_Y_valid == 1] <- (1/M_valid)[Delta_Y_valid == 1] * dnorm(log(M_valid)[Delta_Y_valid == 1], mean = stats::predict(pMX_fit, newdata = data.frame(valid_dat$A, valid_dat$X)[Delta_Y_valid == 1,]), sd = sd(pMX_fit$residuals))
      
      pMXD_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X)[Delta_M_train == 1,])
      # evaluate the fitted density at every non-missing M (Delta_M = 0 participants included), 
      # and set it to 0 outside the observed range of M among Delta_M = 1 in the training data (same rule as haldensify with trim_min = 0)
      range_M_DM <- range(M_train[Delta_M_train == 1])
      it <- index_nomissingM_train
      iv <- index_nomissingM_valid
      pMXDn_A0_train <- rep(NA, n_train)
      pMXDn_A0_train[it] <- (1/M_train[it]) * dnorm(log(M_train[it]), mean = stats::predict(pMXD_fit, newdata = data.frame(A=0, train_dat$X)[it,,drop=FALSE]), sd = sd(pMXD_fit$residuals))
      pMXDn_A0_train[it][M_train[it] < range_M_DM[1] | M_train[it] > range_M_DM[2]] <- 0
      pMXDn_A0_valid <- rep(NA, n_valid)
      pMXDn_A0_valid[iv] <- (1/M_valid[iv]) * dnorm(log(M_valid[iv]), mean = stats::predict(pMXD_fit, newdata = data.frame(A=0, valid_dat$X)[iv,,drop=FALSE]), sd = sd(pMXD_fit$residuals))
      pMXDn_A0_valid[iv][M_valid[iv] < range_M_DM[1] | M_valid[iv] > range_M_DM[2]] <- 0
    }
    
    if(dens_method["pMXZ"] == "hal"){
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
      pMXZn_A0_train <- pMXZn_A1_train <- pMXZn_A_train <- rep(NA, n_train)
      pMXZn_A_train[Delta_Y_train==1] <- stats::predict(pMXZ_fit, new_A = M_train[Delta_Y_train==1], new_W = data.frame(train_dat$A, train_dat$X, train_dat$Z)[Delta_Y_train==1,])
      pMXZn_A0_train[Delta_Y_train==1] <- stats::predict(pMXZ_fit, new_A = M_train[Delta_Y_train==1], new_W = data.frame(A = 0, train_dat$X, train_dat$Z)[Delta_Y_train==1,])
      pMXZn_A1_train[Delta_Y_train==1] <- stats::predict(pMXZ_fit, new_A = M_train[Delta_Y_train==1], new_W = data.frame(A = 1, train_dat$X, train_dat$Z)[Delta_Y_train==1,])
      pMXZn_A0_valid <- pMXZn_A1_valid <- pMXZn_A_valid <- rep(NA, n_valid)
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
      pMXZDn_A_train <- rep(NA, n_train)
      pMXZDn_A_train[index_nomissingM_train] <- stats::predict(pMXZD_fit, new_A = M_train[index_nomissingM_train], new_W = data.frame(train_dat$A, train_dat$X, train_dat$Z)[index_nomissingM_train,], trim_min = 0)
      pMXZDn_A_valid <- rep(NA, n_valid)
      pMXZDn_A_valid[index_nomissingM_valid] <- stats::predict(pMXZD_fit, new_A = M_valid[index_nomissingM_valid], new_W = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z)[index_nomissingM_valid,], trim_min = 0)
    }else if(dens_method["pMXZ"] == "gamlss"){
      # GAMLSS generalized gamma densities p(m|a,x,z,Delta_Y=1) and p(m|a,x,z,Delta_M=1);
      # pMXZD is right-truncated at Delta_M_threshold when it is given
      WZ_train <- data.frame(train_dat$A, train_dat$X, train_dat$Z)
      WZ_valid <- data.frame(valid_dat$A, valid_dat$X, valid_dat$Z)
      dy_t <- Delta_Y_train == 1
      dy_v <- Delta_Y_valid == 1
      it <- index_nomissingM_train
      iv <- index_nomissingM_valid
      pMXZ_fit <- gg_component(M_train[dy_t], WZ_train[dy_t, , drop = FALSE], "pMXZ", gg_opts)
      pMXZn_A0_train <- pMXZn_A1_train <- pMXZn_A_train <- rep(NA, n_train)
      pMXZn_A_train[dy_t] <- gg_density(pMXZ_fit, M_train[dy_t], WZ_train[dy_t, , drop = FALSE], gg_opts)
      pMXZn_A0_train[dy_t] <- gg_density(pMXZ_fit, M_train[dy_t], data.frame(A = 0, train_dat$X, train_dat$Z)[dy_t, , drop = FALSE], gg_opts)
      pMXZn_A1_train[dy_t] <- gg_density(pMXZ_fit, M_train[dy_t], data.frame(A = 1, train_dat$X, train_dat$Z)[dy_t, , drop = FALSE], gg_opts)
      pMXZn_A0_valid <- pMXZn_A1_valid <- pMXZn_A_valid <- rep(NA, n_valid)
      pMXZn_A_valid[dy_v] <- gg_density(pMXZ_fit, M_valid[dy_v], WZ_valid[dy_v, , drop = FALSE], gg_opts)
      pMXZn_A0_valid[dy_v] <- gg_density(pMXZ_fit, M_valid[dy_v], data.frame(A = 0, valid_dat$X, valid_dat$Z)[dy_v, , drop = FALSE], gg_opts)
      pMXZn_A1_valid[dy_v] <- gg_density(pMXZ_fit, M_valid[dy_v], data.frame(A = 1, valid_dat$X, valid_dat$Z)[dy_v, , drop = FALSE], gg_opts)
      
      pMXZD_fit <- gg_component(M_train[Delta_M_train == 1], WZ_train[Delta_M_train == 1, , drop = FALSE], "pMXZD", gg_opts, threshold = Delta_M_threshold)
      pMXZDn_A_train <- rep(NA, n_train)
      pMXZDn_A_train[it] <- gg_density(pMXZD_fit, M_train[it], WZ_train[it, , drop = FALSE], gg_opts, zero_outside = TRUE)
      pMXZDn_A_valid <- rep(NA, n_valid)
      pMXZDn_A_valid[iv] <- gg_density(pMXZD_fit, M_valid[iv], WZ_valid[iv, , drop = FALSE], gg_opts, zero_outside = TRUE)
      density_log <- rbind(density_log, pMXZ_fit$log, pMXZD_fit$log)
    }else{
      pMXZ_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMXZ), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X, train_dat$Z)[Delta_Y_train == 1,])
      pMXZn_A0_train <- pMXZn_A1_train <- pMXZn_A_train <- rep(NA, n_train)
      pMXZn_A_train[Delta_Y_train==1] <- (1/M_train)[Delta_Y_train==1] * dnorm(log(M_train)[Delta_Y_train==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z)[Delta_Y_train==1, ]), sd = sd(pMXZ_fit$residuals))
      pMXZn_A0_train[Delta_Y_train==1] <- (1/M_train)[Delta_Y_train==1] * dnorm(log(M_train)[Delta_Y_train==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 0, train_dat$X, train_dat$Z)[Delta_Y_train==1, ]), sd = sd(pMXZ_fit$residuals))
      pMXZn_A1_train[Delta_Y_train==1] <- (1/M_train)[Delta_Y_train==1] * dnorm(log(M_train)[Delta_Y_train==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 1, train_dat$X, train_dat$Z)[Delta_Y_train==1, ]), sd = sd(pMXZ_fit$residuals))      
      pMXZn_A0_valid <- pMXZn_A1_valid <- pMXZn_A_valid <- rep(NA, n_valid)
      pMXZn_A_valid[Delta_Y_valid==1] <- (1/M_valid)[Delta_Y_valid==1] * dnorm(log(M_valid)[Delta_Y_valid==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z)[Delta_Y_valid==1, ]), sd = sd(pMXZ_fit$residuals))
      pMXZn_A0_valid[Delta_Y_valid==1] <- (1/M_valid)[Delta_Y_valid==1] * dnorm(log(M_valid)[Delta_Y_valid==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z)[Delta_Y_valid==1, ]), sd = sd(pMXZ_fit$residuals))
      pMXZn_A1_valid[Delta_Y_valid==1] <- (1/M_valid)[Delta_Y_valid==1] * dnorm(log(M_valid)[Delta_Y_valid==1], mean = stats::predict(pMXZ_fit, newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z)[Delta_Y_valid==1, ]), sd = sd(pMXZ_fit$residuals))      
      
      pMXZD_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMXZ), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X, train_dat$Z)[Delta_M_train == 1,])
      # evaluate the fitted density at every non-missing M; 0 outside the observed range of M among Delta_M = 1 in the training data
      range_M_DM <- range(M_train[Delta_M_train == 1])
      it <- index_nomissingM_train
      iv <- index_nomissingM_valid
      pMXZDn_A_train <- rep(NA, n_train)
      pMXZDn_A_train[it] <- (1/M_train[it]) * dnorm(log(M_train[it]), mean = stats::predict(pMXZD_fit, newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z)[it,,drop=FALSE]), sd = sd(pMXZD_fit$residuals))
      pMXZDn_A_train[it][M_train[it] < range_M_DM[1] | M_train[it] > range_M_DM[2]] <- 0
      pMXZDn_A_valid <- rep(NA, n_valid)
      pMXZDn_A_valid[iv] <- (1/M_valid[iv]) * dnorm(log(M_valid[iv]), mean = stats::predict(pMXZD_fit, newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z)[iv,,drop=FALSE]), sd = sd(pMXZD_fit$residuals))
      pMXZDn_A_valid[iv][M_valid[iv] < range_M_DM[1] | M_valid[iv] > range_M_DM[2]] <- 0
    }
  }
  
  # density-ratio weights pMXD / pMXZ that enter the efficient influence function (validation fold)
  if(M_indicator){
    A_valid <- valid_dat$A$A
    weight_diagnostics <- list(
      A0 = (pMXDn_A0_valid / pMXZn_A0_valid)[Delta_Y_valid == 1 & A_valid == 0],
      A1 = (pMXDn_A0_valid / pMXZn_A1_valid)[Delta_Y_valid == 1 & A_valid == 1]
    )
  }
  
  # define parameters for storing results
  est <- matrix(nrow = 2, ncol = p)
  eif_mat <- vector("list", p)
  for(j in 1:p){
    # fit outcome regression
    mu_MXZn_A0_train <- mu_MXZn_A1_train <- mu_MXZn_A_train <- rep(NA, n_train)
    mu_MXZn_A0_valid <- mu_MXZn_A1_valid <- mu_MXZn_A_valid <- rep(NA, n_valid)
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
                                                SL.library = SL_mu_AMXZ,
                                                method = tmp_method.CC_LS(),
                                                control = list(saveCVFitLibrary = TRUE))
      mu_MXZn_A0_train[index_nomissingM_train] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, train_dat$M, train_dat$X, train_dat$Z)[index_nomissingM_train,], onlySL = TRUE)[[1]]
      mu_MXZn_A1_train[index_nomissingM_train] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, train_dat$M, train_dat$X, train_dat$Z)[index_nomissingM_train,], onlySL = TRUE)[[1]]
      mu_MXZn_A0_valid[index_nomissingM_valid] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, valid_dat$M, valid_dat$X, valid_dat$Z)[index_nomissingM_valid,], onlySL = TRUE)[[1]]
      mu_MXZn_A1_valid[index_nomissingM_valid] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, valid_dat$M, valid_dat$X, valid_dat$Z)[index_nomissingM_valid,], onlySL = TRUE)[[1]]
      mu_MXZn_A_train[index_nomissingM_train] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(train_dat$A, train_dat$M, train_dat$X, train_dat$Z)[index_nomissingM_train,], onlySL = TRUE)[[1]]
      mu_MXZn_A_valid[index_nomissingM_valid] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(valid_dat$A, valid_dat$M, valid_dat$X, valid_dat$Z)[index_nomissingM_valid,], onlySL = TRUE)[[1]]
    }
    
    if(M_indicator){
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
                                                  SL.library = SL_eta_AXZ,
                                                  method = tmp_method.CC_LS(),
                                                  control = list(saveCVFitLibrary = TRUE))
        eta_AXZn_A_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z), onlySL = TRUE)[[1]]
        eta_AXZn_A0_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, train_dat$X, train_dat$Z), onlySL = TRUE)[[1]]
        eta_AXZn_A1_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, train_dat$X, train_dat$Z), onlySL = TRUE)[[1]]
        
        eta_AXZn_A_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z), onlySL = TRUE)[[1]]
        eta_AXZn_A0_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z), onlySL = TRUE)[[1]]
        eta_AXZn_A1_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z), onlySL = TRUE)[[1]]
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
      }else if(!is.null(SL_eta_AXM)){
        set.seed(seed)
        eta_AXMn_A0_train <- eta_AXMn_A1_train <- rep(NA, n_train)
        eta_AXMn_A0_valid <- eta_AXMn_A1_valid <- rep(NA, n_valid)
        # since observations with Delta=0 will not be considered in the next step
        # so add restrictions here to better depict data with Delta=1
        eta_AXM_fit <- SuperLearner::SuperLearner(Y = mu_pseudo_A_star_train[Delta_Y_train==1], 
                                                  X = data.frame(train_dat$A, train_dat$M, train_dat$X)[Delta_Y_train==1,],
                                                  family = gaussian(),
                                                  SL.library = SL_eta_AXM,
                                                  method = tmp_method.CC_LS(),
                                                  control = list(saveCVFitLibrary = TRUE))
        eta_AXMn_A0_train[index_nomissingM_train] <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, train_dat$M, train_dat$X)[index_nomissingM_train,], onlySL = TRUE)[[1]]
        eta_AXMn_A1_train[index_nomissingM_train] <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, train_dat$M, train_dat$X)[index_nomissingM_train,], onlySL = TRUE)[[1]]
        
        eta_AXMn_A0_valid[index_nomissingM_valid] <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, valid_dat$M, valid_dat$X)[index_nomissingM_valid,], onlySL = TRUE)[[1]]
        eta_AXMn_A1_valid[index_nomissingM_valid] <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, valid_dat$M, valid_dat$X)[index_nomissingM_valid,], onlySL = TRUE)[[1]]
      }else{
        if(!requireNamespace("hal9001", quietly = TRUE)){
          stop("Package 'hal9001' is required when no SuperLearner library or glm formula is given for eta_AXM.")
        }
        set.seed(seed)
        eta_AXMn_A0_train <- eta_AXMn_A1_train <- rep(NA, n_train)
        eta_AXMn_A0_valid <- eta_AXMn_A1_valid <- rep(NA, n_valid)
        eta_AXM_fit <- hal9001::fit_hal(X = data.frame(train_dat$A, train_dat$M, train_dat$X)[Delta_Y_train==1,], Y = mu_pseudo_A_star_train[Delta_Y_train==1]) 
        eta_AXMn_A0_train[index_nomissingM_train] <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 0, train_dat$M, train_dat$X)[index_nomissingM_train,,drop=FALSE])
        eta_AXMn_A1_train[index_nomissingM_train] <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 1, train_dat$M, train_dat$X)[index_nomissingM_train,,drop=FALSE])
        
        eta_AXMn_A0_valid[index_nomissingM_valid] <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 0, valid_dat$M, valid_dat$X)[index_nomissingM_valid,,drop=FALSE])
        eta_AXMn_A1_valid[index_nomissingM_valid] <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 1, valid_dat$M, valid_dat$X)[index_nomissingM_valid,,drop=FALSE])
      }
      
    }else{
      eta_AXZn_A_train = mu_MXZn_A_train
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
                                           SL.library = SL_xi_AX,
                                           method = tmp_method.CC_LS(),
                                           control = list(saveCVFitLibrary = TRUE))
      xi_AXn_A0 <- stats::predict(xi_fit, newdata = data.frame(A = 0, valid_dat$X), onlySL = TRUE)[[1]]
      xi_AXn_A1 <- stats::predict(xi_fit, newdata = data.frame(A = 1, valid_dat$X), onlySL = TRUE)[[1]]
    }
    
    # calculate efficient influence function
    if(M_indicator){
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
    }else{
      eif_A0 <- make_full_data_eif_easy(a = 0, A = valid_dat$A$A, Delta_Y = Delta_Y_valid, Y = valid_dat$Y[,j], 
                                        gA = 1 - gAn_1, gDY_AXZ = gDYn_1_A0XZ_valid,
                                        xi_AX = xi_AXn_A0, mu = mu_MXZn_A0_valid)
      
      eif_A1 <- make_full_data_eif_easy(a = 1, A = valid_dat$A$A, Delta_Y = Delta_Y_valid, Y = valid_dat$Y[,j], 
                                        gA = gAn_1, gDY_AXZ = gDYn_1_A1XZ_valid,
                                        xi_AX = xi_AXn_A1, mu = mu_MXZn_A1_valid)
    }
    
    # one-step estimator
    est[,j] <- c(mean(xi_AXn_A0) + mean(eif_A0), mean(xi_AXn_A1) + mean(eif_A1))
    
    # covariance matrix
    eif_mat[[j]] <- cbind(eif_A0, eif_A1)
    colnames(eif_mat[[j]]) <- c("eif_A0", "eif_A1")
  }

  # output
  out <- list(est = est,
              eif_mat = eif_mat,
              density_log = density_log,
              weights = weight_diagnostics
  )
  
  return(out)
}

#' Non-parametric efficient motion-controlled brain phenotype estimators
#' 
#' @details Compute non-parametric efficient motion-controlled brain phenotype estimators (MoCo) using one-step cross-fitting
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
#'                    - \code{gD}: GLM formula for estimating the probability P(Delta_M = 1 | A, X).
#'                    - \code{mu_AMXZ}: GLM formula for estimating the outcome regression E(Y | A, M, X, Z).
#'                    - \code{xi_AX}: GLM formula for estimating E(eta_AXZ | A=a, X).
#'                    - \code{eta_AXM}: GLM formula for estimating E(Q * p / r | A=a, M, X).
#'                    - \code{pMX}: GLM formula for estimating p(m | a, x) and p(m | a, x, Delta_M = 1), assuming M follows a normal distribution.
#'                    - \code{pMXZ}: GLM formula for estimating p(m | a, x, z) and p(m | a, x, z, Delta_M = 1), assuming M follows a normal distribution.
#'                    
#' @param HAL_pMX Specifies whether to estimate p(m | a, x, Delta_Y = 1) and p(m | a, x, Delta_M=1) using the highly adaptive lasso conditional density estimation method. 
#' Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMX} option in \code{glm_formula}, such as \code{pMX = "."}.
#' @param HAL_pMXZ Specifies whether to estimate p(m | a, x, z, Delta_Y = 1) and p(m | a, x, z, Delta_M=1) using the highly adaptive lasso conditional density estimation method. 
#' Defaults to \code{TRUE}. If set to \code{FALSE}, please specify the \code{pMXZ} option in \code{glm_formula}, such as \code{pMXZ = "."}.
#'  
#' @param HAL_options Additional options for highly adaptive lasso (HAL) method.
#'                   - \code{max_degree}: The highest order of interaction terms for generating basis functions (passed to \code{haldensify}).
#'                   - \code{lambda_seq}: A numeric sequence of values for the regularization parameter of Lasso regression (passed to \code{haldensify}).
#'                   - \code{num_knots}: The maximum number of knot points (i.e., bins) for any covariate for generating basis functions (passed to \code{haldensify}).
#' 
#' @param cv_folds a numeric value indicates the number(s) of partitions in cross-fitting. Defaults to 5. 
#' 
#' @param density_method,Delta_M_threshold,gamlss_options Motion density method, the threshold that defines \code{Delta_M}, and GAMLSS options. See \code{\link{moco}}.
#' @param seed Random seed set before each SuperLearner fit (one value of \code{seed_rng} in \code{moco()}).
#' @param ... Not used.
#' 
#' @import SuperLearner 
#' @import haldensify 
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
      gDY_X = NULL,
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
    HAL_options = list(max_degree = 3, lambda_seq = exp(seq(-1, -10, length = 100)), num_knots = c(1000, 500, 250)),
    density_method = NULL,
    Delta_M_threshold = NULL,
    gamlss_options = NULL,
    seed = 1, 
    cv_folds = 5,
    ...
){
  # number of participants
  n <- nrow(Y)
  # number of edges
  p <- ncol(Y)
  # whether motion is available
  M_indicator <- (!is.null(M))
  
  if(!M_indicator){
    # missingness indicator as in Nebel
    Delta_Y <- ifelse(Delta_M + Delta_Y == 2, 1, 0)
    # impute M to be all 1's for computation convenience
    M <- rep(1, n)
    # M[Delta_Y == 0] <- NA
  }
  
  # if X or Z is null
  if(is.null(X)){X <- rep(1, n)}
  if(is.null(Z)){Z <- rep(1, n)}
  
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
  eif_mat <- vector("list", p)
  for(j in 1:p){
    eif_mat[[j]] <- matrix(nrow = n, ncol = 2)
  }
  
  # density fit log and density-ratio weights, collected over folds
  density_log <- NULL
  weights_A0 <- weights_A1 <- NULL
  
  # cross fitting
  for(i in 1:cv_folds){
    rslt_cross <- fit_mechanism(
      train_dat = lapply(dat, function(x){x[-folds[[i]], , drop = FALSE]}), 
      valid_dat = lapply(dat, function(x){x[folds[[i]], , drop = FALSE]}),
      M_indicator = M_indicator,
      SL_library = SL_library,
      SL_library_customize = SL_library_customize,
      glm_formula = glm_formula,
      HAL_pMX = HAL_pMX,
      HAL_pMXZ = HAL_pMXZ,
      HAL_options = HAL_options,
      density_method = density_method,
      Delta_M_threshold = Delta_M_threshold,
      gamlss_options = gamlss_options,
      seed = seed
    )

    # store results
    est <- est + rslt_cross$est
    if(!is.null(rslt_cross$density_log)){
      density_log <- rbind(density_log, data.frame(fold = i, rslt_cross$density_log))
    }
    if(!is.null(rslt_cross$weights)){
      weights_A0 <- c(weights_A0, rslt_cross$weights$A0)
      weights_A1 <- c(weights_A1, rslt_cross$weights$A1)
    }
    for(j in 1:p){
      eif_mat[[j]][folds[[i]],] <- rslt_cross$eif_mat[[j]]
    }
  }
  
  # calculate one-step estimator
  est <- est/cv_folds 
  
  # motion-adjusted associations
  adj_association <- est[2, ] - est[1, ]
  
  # calculate covariance matrix
  cov_mat <- vector("list", p)
  for(j in 1:p){
    colnames(eif_mat[[j]]) <- c('est_A0', 'est_A1')
    cov_mat[[j]] <- cov(eif_mat[[j]])/n
  }
  
  # summary of the density-ratio weights pMXD / pMXZ in each arm
  weight_diagnostics <- NULL
  if(M_indicator){
    weight_diagnostics <- rbind(A0 = weight_summary(weights_A0), A1 = weight_summary(weights_A1))
  }
  
  # output
  out <- list(est = est,
              adj_association = adj_association,
              eif_mat = eif_mat,
              cov_mat = cov_mat,
              density_log = density_log,
              weight_diagnostics = weight_diagnostics
  )
}

