# Internal helpers for the generalized-Gamma motion-density implementation.
#
# The default specification uses penalized B-splines for continuous
# covariates in both the location and scale submodels, includes categorical
# covariates linearly, includes the exposure/group indicator in the location
# submodel, and fixes the shape parameter. These helpers are intentionally
# dataset agnostic; study-specific variable names belong in analysis scripts.

.moco_bt <- function(x) {
  paste0("`", gsub("`", "\\`", x, fixed = TRUE), "`")
}

.moco_join_terms <- function(x) {
  if (length(x) == 0L) "1" else paste(x, collapse = " + ")
}

.moco_default_pb_terms <- function(vars) {
  if (length(vars) == 0L) return(character(0))
  vapply(
    vars,
    function(v) paste0("gamlss::pb(", .moco_bt(v), ")"),
    character(1)
  )
}

.moco_linear_terms <- function(vars) {
  if (length(vars) == 0L) return(character(0))
  vapply(vars, .moco_bt, character(1))
}

.moco_make_default_gamlss_formulas <- function(
    x_names,
    z_names = character(0),
    continuous_X = character(0),
    continuous_Z = character(0)
) {
  all_names <- c(x_names, z_names)
  continuous <- unique(c(continuous_X, continuous_Z))
  continuous <- intersect(continuous, all_names)
  discrete <- setdiff(all_names, continuous)

  pb_cont <- .moco_default_pb_terms(continuous)
  lin_disc <- .moco_linear_terms(discrete)

  mu_terms <- c("A", pb_cont, lin_disc)
  sigma_terms <- c(pb_cont, lin_disc)

  list(
    mu = stats::as.formula(
      paste("M ~", .moco_join_terms(mu_terms))
    ),
    sigma = stats::as.formula(
      paste("~", .moco_join_terms(sigma_terms))
    ),
    nu = ~ 1
  )
}

.moco_match_gamlss_optimizer <- function(gamlss_optimizer) {
  if (length(gamlss_optimizer) != 1L || is.na(gamlss_optimizer)) {
    stop("gamlss_optimizer must be one of: 'RS', 'mixed', or 'CG'.")
  }

  key <- tolower(as.character(gamlss_optimizer))
  matched <- unname(c(rs = "RS", mixed = "mixed", cg = "CG")[key])

  if (length(matched) != 1L || is.na(matched)) {
    stop("gamlss_optimizer must be one of: 'RS', 'mixed', or 'CG'.")
  }

  matched
}

.moco_gamlss_optimizer_label <- function(gamlss_optimizer) {
  switch(
    .moco_match_gamlss_optimizer(gamlss_optimizer),
    RS = "RS()",
    mixed = "mixed(1,50)",
    CG = "CG()"
  )
}

.moco_fit_gamlss <- function(
    formula,
    sigma.formula,
    nu.formula,
    data,
    control,
    gamlss_optimizer = "RS",
    model_label = "GAMLSS"
) {
  gamlss_optimizer <- .moco_match_gamlss_optimizer(gamlss_optimizer)

  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (!"M" %in% names(data)) stop(model_label, ": fitting data do not contain M.")
  if (nrow(data) == 0L) stop(model_label, ": no observations remain for fitting.")
  if (anyNA(data)) stop(model_label, ": fitting data contain missing values.")
  if (!is.numeric(data$M)) stop(model_label, ": M must be numeric for the GG family.")

  numeric_columns <- vapply(data, is.numeric, logical(1))
  if (any(vapply(data[numeric_columns], function(x) any(!is.finite(x)), logical(1)))) {
    stop(model_label, ": fitting data contain non-finite numeric values.")
  }
  if (any(data$M <= 0)) {
    stop(model_label, ": the GG family requires strictly positive M values.")
  }

  fit <- tryCatch(
    switch(
      gamlss_optimizer,
      RS = gamlss::gamlss(
        formula = formula,
        sigma.formula = sigma.formula,
        nu.formula = nu.formula,
        family = "GG",
        data = data,
        method = RS(),
        control = control
      ),
      mixed = gamlss::gamlss(
        formula = formula,
        sigma.formula = sigma.formula,
        nu.formula = nu.formula,
        family = "GG",
        data = data,
        method = mixed(1, 50),
        control = control
      ),
      CG = gamlss::gamlss(
        formula = formula,
        sigma.formula = sigma.formula,
        nu.formula = nu.formula,
        family = "GG",
        data = data,
        method = CG(),
        control = control
      )
    ),
    error = function(e) {
      stop(
        model_label, " failed with optimizer ",
        .moco_gamlss_optimizer_label(gamlss_optimizer), ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (!isTRUE(fit$converged)) {
    stop(
      model_label, " did not converge with optimizer ",
      .moco_gamlss_optimizer_label(gamlss_optimizer), ".",
      call. = FALSE
    )
  }
  if (length(fit$G.deviance) != 1L || !is.finite(fit$G.deviance)) {
    stop(model_label, ": fitted global deviance is not finite.", call. = FALSE)
  }

  bad_parameter <- vapply(
    fit$parameters,
    function(parameter) {
      fitted_values <- fit[[paste0(parameter, ".fv")]]
      is.null(fitted_values) || any(!is.finite(fitted_values))
    },
    logical(1)
  )
  if (any(bad_parameter)) {
    stop(
      model_label, ": non-finite fitted values for parameter(s): ",
      paste(fit$parameters[bad_parameter], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  attr(fit, "moco_gamlss_optimizer") <- .moco_gamlss_optimizer_label(gamlss_optimizer)
  fit
}

##########################################################
# one_step function with gamlss motion density
# similarly, we can use gamlss in one_step_cross function
##########################################################

one_step_gamlss <- function(
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
    GAMLSS_pMX = FALSE,
    GAMLSS_pMXZ = FALSE,
    gamlss_formula = list(
      pMX_mu = NULL,
      pMX_sigma = ~ 1,
      pMX_nu = ~ 1,
      pMXZ_mu = NULL,
      pMXZ_sigma = ~ 1,
      pMXZ_nu = ~ 1
    ),
    gamlss_family = "GG",
    gamlss_optimizer = "RS",
    GAMLSS_BIC_select = FALSE,
    gamlss_continuous_X = character(0),
    gamlss_continuous_Z = character(0),
    gamlss_bic_candidates = c("linear", "pb_mu", "pb_mu_sigma"),
    gamlss_bic_trace = FALSE,
    n.cyc = 300,
    HAL_pMX = FALSE,
    HAL_pMXZ = FALSE,
    HAL_options = list(
      max_degree = 3,
      lambda_seq = exp(seq(-1, -10, length = 100)),
      num_knots = c(1000, 500, 250)
    ),
    seed = 1,
    ...
){
  if (!identical(gamlss_family, "GG")) {
    stop("The current GAMLSS density implementation supports gamlss_family = \"GG\" only.")
  }
  if (isTRUE(GAMLSS_pMX) && isTRUE(HAL_pMX)) {
    stop("Choose only one pMX density method: GAMLSS or HAL.")
  }
  if (isTRUE(GAMLSS_pMXZ) && isTRUE(HAL_pMXZ)) {
    stop("Choose only one pMXZ density method: GAMLSS or HAL.")
  }
  # number of participants
  n <- nrow(Y)
  # number of outcomes
  p <- ncol(Y)
  
  if(!is.null(thresh)){
    # if thresh is not null, collapse M into dummy variables based on truncated level thresh
    Delta_M <- as.numeric(M < thresh)
    Delta_M[is.na(M)] <- 0
  }
  
  # change X and Z to a dataframe if is a vector
  if(is.null(dim(X))){X <- data.frame(X = X)}
  if(is.null(dim(Z))){Z <- data.frame(Z = Z)}
  
  # Store BIC model-selection diagnostics for non-cross-fit analyses.
  gamlss_selection <- list()
  
  # Explicitly named data frames prevent formula variables/splines from
  # being resolved from an outer environment.
  A_num <- as.numeric(A)
  pMX_all_data <- data.frame(M = M, A = A_num, X, check.names = FALSE)
  pMXZ_all_data <- data.frame(M = M, A = A_num, X, Z, check.names = FALSE)
  
  # GAMLSS parameter prediction through the standard S3 predict generic.
  # This uses stats::predict() for mu, sigma, and nu.
  predict_gamlss_parameters_nocv <- function(fit, newdata, fit_data) {
    list(
      mu = as.numeric(stats::predict(
        fit, what = "mu", type = "response",
        newdata = newdata, data = fit_data
      )),
      sigma = as.numeric(stats::predict(
        fit, what = "sigma", type = "response",
        newdata = newdata, data = fit_data
      )),
      nu = as.numeric(stats::predict(
        fit, what = "nu", type = "response",
        newdata = newdata, data = fit_data
      ))
    )
  }
  
  # non-missing motion index
  index_nomissingM <- which(!is.na(M))
  
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
  
  ####################################
  # fit regression for propensity score
  # exposure or group indicator, conditional on baseline covariates
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
  
  ####################################
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
  
  ############################################
  # fit missingness indicator regression
  if(sum(Delta_Y) == n){
    gDYn_1_AX <- rep(1, n)
    gDYn_1_AXZ <- gDYn_1_A0XZ <- gDYn_1_A1XZ <- rep(1, n)
  }else{
    # probability of non-missing conditioning on exposure or group indicator A, baseline covariates X and post-exposure covariates Z
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
    
    # probability of non-missing conditioning on exposure or group indicator A, baseline covariates X and post-exposure covariates Z
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
  
  ########################################
  # Motion density
  if(HAL_pMX){
    # density estimation for mediator, conditioning on exposure or group indicator and baseline covariates p(m|a,x)
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
    
    # density estimation for mediator, conditioning on exposure or group indicator and baseline covariates p(m|0,x,Delta_M=1)
    pMXD_fit <- haldensify::haldensify(
      A = M[Delta_M == 1],
      W = data.frame(A, X)[Delta_M == 1,],
      max_degree = HAL_options$max_degree,
      lambda_seq = HAL_options$lambda_seq,
      num_knots = HAL_options$num_knots
    )
    pMXDn_A0 <- rep(NA, n)
    pMXDn_A0[index_nomissingM] <- stats::predict(pMXD_fit, new_A = M[index_nomissingM], new_W = data.frame(A = 0, X)[index_nomissingM,], trim_min = 0)
  }else if(GAMLSS_pMX){
    ctrl <- gamlss::gamlss.control(n.cyc = n.cyc, c.crit = 1e-6, trace = isTRUE(gamlss_bic_trace))
    
    # --------------------------------------------------------------
    # p(M | A, X)
    #
    # Fit on all observations with positive observed M, matching the
    # original MoCo density definition. Density values needed by the
    # outcome-missingness part are evaluated where Delta_Y == 1.
    # --------------------------------------------------------------
    idx_all_M <- which(!is.na(M) & is.finite(M) & M > 0)
    idx_Y <- which(Delta_Y == 1 & !is.na(M) & is.finite(M) & M > 0)
    pMX_fit_data <- pMX_all_data[idx_all_M, , drop = FALSE]
    if (nrow(pMX_fit_data) == 0L) stop("No observations remain for pMX GAMLSS fitting.")

    pMX_spec <- .moco_make_default_gamlss_formulas(
      x_names = names(X),
      continuous_X = gamlss_continuous_X
    )
    pMX_fit <- .moco_fit_gamlss(
      formula = pMX_spec$mu,
      sigma.formula = pMX_spec$sigma,
      nu.formula = pMX_spec$nu,
      data = pMX_fit_data,
      control = ctrl,
      gamlss_optimizer = gamlss_optimizer,
      model_label = "pMX GAMLSS"
    )
    gamlss_selection$pMX <- list(
      method = "default_pb_mu_sigma_nu1",
      optimizer = .moco_gamlss_optimizer_label(gamlss_optimizer),
      formulas = pMX_spec
    )

    pMXn_A <- rep(NA_real_, n)
    pMX_predict <- predict_gamlss_parameters_nocv(
      pMX_fit, pMX_all_data[idx_Y, , drop = FALSE], pMX_fit_data
    )
    pMXn_A[idx_Y] <- gamlss.dist::dGG(
      M[idx_Y], mu = pMX_predict$mu, sigma = pMX_predict$sigma, nu = pMX_predict$nu
    )
    
    # --------------------------------------------------------------
    # p(M | A = 0, X, Delta_M = 1)
    # --------------------------------------------------------------
    idx_M <- which(Delta_M == 1 & !is.na(M) & is.finite(M) & M > 0)
    pMXD_fit_data <- pMX_all_data[idx_M, , drop = FALSE]
    if (nrow(pMXD_fit_data) == 0L) stop("No observations remain for pMXD GAMLSS fitting.")

    pMXD_spec <- .moco_make_default_gamlss_formulas(
      x_names = names(X),
      continuous_X = gamlss_continuous_X
    )
    pMXD_fit <- .moco_fit_gamlss(
      formula = pMXD_spec$mu,
      sigma.formula = pMXD_spec$sigma,
      nu.formula = pMXD_spec$nu,
      data = pMXD_fit_data,
      control = ctrl,
      gamlss_optimizer = gamlss_optimizer,
      model_label = "pMXD GAMLSS"
    )
    gamlss_selection$pMXD <- list(
      method = "default_pb_mu_sigma_nu1",
      optimizer = .moco_gamlss_optimizer_label(gamlss_optimizer),
      formulas = pMXD_spec
    )

    pMXDn_A0 <- rep(NA_real_, n)
    pMXDn_A0[Delta_M == 0] <- 0
    pMXD_new <- pMX_all_data[idx_M, , drop = FALSE]
    pMXD_new$A <- 0
    pMXD_predict <- predict_gamlss_parameters_nocv(pMXD_fit, pMXD_new, pMXD_fit_data)
    pMXDn_A0[idx_M] <- gamlss.dist::dGG(
      M[idx_M], mu = pMXD_predict$mu, sigma = pMXD_predict$sigma, nu = pMXD_predict$nu
    )
    
  }else{
    pMX_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M), A, X)[Delta_Y == 1,])
    pMXn_A <- rep(NA, n)
    pMXn_A[Delta_Y == 1] <- (1/M)[Delta_Y == 1] * dnorm(log(M)[Delta_Y == 1], mean = stats::predict(pMX_fit, newdata = data.frame(A, X)[Delta_Y == 1,]), sd = sd(pMX_fit$residuals))
    
    pMXD_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M), A, X)[Delta_M == 1,])
    pMXDn_A0 <- rep(NA, n)
    pMXDn_A0[Delta_M==0] <- 0
    pMXDn_A0[Delta_M==1] <- (1/M)[Delta_M == 1] * dnorm(log(M)[Delta_M == 1], mean = stats::predict(pMXD_fit, newdata = data.frame(A=0,X)[Delta_M == 1,]), sd = sd(pMXD_fit$residuals))
  }
  
  if(HAL_pMXZ){
    # density estimation for mediator, conditioning on exposure or group indicator, binary mediator,
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
    
    # density estimation for mediator, conditioning on exposure or group indicator, binary mediator,
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
    pMXZDn_A[index_nomissingM] <- stats::predict(pMXZD_fit, new_A = M[index_nomissingM], new_W = data.frame(A, X, Z)[index_nomissingM,], trim_min = 0)
  } else if(GAMLSS_pMXZ){
    ctrl <- gamlss::gamlss.control(n.cyc = n.cyc, c.crit = 1e-6, trace = isTRUE(gamlss_bic_trace))
    
    # --------------------------------------------------------------
    # p(M | A, X, Z)
    #
    # Fit on all observations with positive observed M.
    # --------------------------------------------------------------
    idx_all_M <- which(!is.na(M) & is.finite(M) & M > 0)
    idx_Y <- which(Delta_Y == 1 & !is.na(M) & is.finite(M) & M > 0)
    pMXZ_fit_data <- pMXZ_all_data[idx_all_M, , drop = FALSE]
    if (nrow(pMXZ_fit_data) == 0L) stop("No observations remain for pMXZ GAMLSS fitting.")

    pMXZ_spec <- .moco_make_default_gamlss_formulas(
      x_names = names(X),
      z_names = names(Z),
      continuous_X = gamlss_continuous_X,
      continuous_Z = gamlss_continuous_Z
    )
    pMXZ_fit <- .moco_fit_gamlss(
      formula = pMXZ_spec$mu,
      sigma.formula = pMXZ_spec$sigma,
      nu.formula = pMXZ_spec$nu,
      data = pMXZ_fit_data,
      control = ctrl,
      gamlss_optimizer = gamlss_optimizer,
      model_label = "pMXZ GAMLSS"
    )
    gamlss_selection$pMXZ <- list(
      method = "default_pb_mu_sigma_nu1",
      optimizer = .moco_gamlss_optimizer_label(gamlss_optimizer),
      formulas = pMXZ_spec
    )

    pMXZn_A <- pMXZn_A0 <- pMXZn_A1 <- rep(NA_real_, n)
    new_obs <- pMXZ_all_data[idx_Y, , drop = FALSE]
    pred_A <- predict_gamlss_parameters_nocv(pMXZ_fit, new_obs, pMXZ_fit_data)
    new_A0 <- new_obs; new_A0$A <- 0
    new_A1 <- new_obs; new_A1$A <- 1
    pred_A0 <- predict_gamlss_parameters_nocv(pMXZ_fit, new_A0, pMXZ_fit_data)
    pred_A1 <- predict_gamlss_parameters_nocv(pMXZ_fit, new_A1, pMXZ_fit_data)
    
    pMXZn_A[idx_Y] <- gamlss.dist::dGG(M[idx_Y], mu = pred_A$mu, sigma = pred_A$sigma, nu = pred_A$nu)
    pMXZn_A0[idx_Y] <- gamlss.dist::dGG(M[idx_Y], mu = pred_A0$mu, sigma = pred_A0$sigma, nu = pred_A0$nu)
    pMXZn_A1[idx_Y] <- gamlss.dist::dGG(M[idx_Y], mu = pred_A1$mu, sigma = pred_A1$sigma, nu = pred_A1$nu)
    
    # --------------------------------------------------------------
    # p(M | A, X, Z, Delta_M = 1)
    # --------------------------------------------------------------
    idx_M <- which(Delta_M == 1 & !is.na(M) & is.finite(M) & M > 0)
    pMXZD_fit_data <- pMXZ_all_data[idx_M, , drop = FALSE]
    if (nrow(pMXZD_fit_data) == 0L) stop("No observations remain for pMXZD GAMLSS fitting.")

    pMXZD_spec <- .moco_make_default_gamlss_formulas(
      x_names = names(X),
      z_names = names(Z),
      continuous_X = gamlss_continuous_X,
      continuous_Z = gamlss_continuous_Z
    )
    pMXZD_fit <- .moco_fit_gamlss(
      formula = pMXZD_spec$mu,
      sigma.formula = pMXZD_spec$sigma,
      nu.formula = pMXZD_spec$nu,
      data = pMXZD_fit_data,
      control = ctrl,
      gamlss_optimizer = gamlss_optimizer,
      model_label = "pMXZD GAMLSS"
    )
    gamlss_selection$pMXZD <- list(
      method = "default_pb_mu_sigma_nu1",
      optimizer = .moco_gamlss_optimizer_label(gamlss_optimizer),
      formulas = pMXZD_spec
    )

    pMXZDn_A <- rep(NA_real_, n)
    pMXZDn_A[Delta_M == 0] <- 0
    pred_M <- predict_gamlss_parameters_nocv(
      pMXZD_fit, pMXZ_all_data[idx_M, , drop = FALSE], pMXZD_fit_data
    )
    pMXZDn_A[idx_M] <- gamlss.dist::dGG(
      M[idx_M], mu = pred_M$mu, sigma = pred_M$sigma, nu = pred_M$nu
    )
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
    mu_MXZn_A0 <- mu_MXZn_A1 <- mu_MXZn_A <- rep(NA, n)
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
      mu_MXZn_A0[index_nomissingM] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, M, X, Z)[index_nomissingM,])[[1]]
      mu_MXZn_A1[index_nomissingM] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, M, X, Z)[index_nomissingM,])[[1]]
      mu_MXZn_A[index_nomissingM] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A, M, X, Z)[index_nomissingM,])[[1]]
    }
    
    # fit additional pseudo-outcome regression for mu_MXZ*pMXD/pMXZD, conditioning on exposure or group indicator, binary mediator and baseline covariates
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
    
    # fit pseudo-outcome regression for Qd, conditioning on exposure or group indicator, mediator and baseline covariates
    mu_pseudo_A_star <- mu_MXZn_A*pMXn_A/pMXZn_A*gDYn_1_AX/gDYn_1_AXZ #warning here
    # fit the regression
    if(!is.null(glm_formula$eta_AXM)){
      eta_AXM_fit <- stats::glm(paste0("mu_pseudo_A_star ~ ", glm_formula$eta_AXM), family = gaussian(),
                                data = data.frame(mu_pseudo_A_star = mu_pseudo_A_star, A, M, X)[Delta_Y == 1,])
      eta_AXMn_A0 <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, M, X))
      eta_AXMn_A1 <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, M, X))
    }else if(!is.null(SL_eta_AXM)){
      set.seed(seed)
      eta_AXMn_A0 <- eta_AXMn_A1 <- rep(NA, n)
      # since observations with Delta_M=0 will not be considered in the next step
      # so add restrictions here to better depict data with Delta_M=1 # [Delta_M==1]
      eta_AXM_fit <- SuperLearner::SuperLearner(Y = mu_pseudo_A_star[Delta_Y == 1], X = data.frame(A, M, X)[Delta_Y == 1,],
                                                family = gaussian(),
                                                SL.library = SL_eta_AXM,
                                                method = tmp_method.CC_LS(),
                                                control = list(saveCVFitLibrary = TRUE))
      eta_AXMn_A0[index_nomissingM] <- stats::predict(eta_AXM_fit, type = "response", newdata = data.frame(A = 0, M, X)[index_nomissingM,])[[1]]
      eta_AXMn_A1[index_nomissingM] <- stats::predict(eta_AXM_fit, type = "response", newdata = data.frame(A = 1, M, X)[index_nomissingM,])[[1]]
    }else{
      set.seed(seed)
      eta_AXM_fit <- hal9001::fit_hal(X = data.frame(A, M, X)[Delta_Y == 1,], Y = mu_pseudo_A_star[Delta_Y == 1])
      eta_AXMn_A0 <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 0, M, X))
      eta_AXMn_A1 <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 1, M, X))
    }
    
    # fit pseudo-outcome regression for eta_AXZ, conditioning on exposure or group indicator and baseline covariates
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
              cov_mat = cov_mat,
              gamlss_selection = gamlss_selection
  )
  
  return(out)
}

#######################################
# one step cross with gamlss
#######################################

fit_mechanism_gamlss <- function(
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
    GAMLSS_pMX = FALSE,
    GAMLSS_pMXZ = FALSE,
    gamlss_formula = list(
      pMX_mu = NULL,
      pMX_sigma = ~ 1,
      pMX_nu = ~ 1,
      pMXZ_mu = NULL,
      pMXZ_sigma = ~ 1,
      pMXZ_nu = ~ 1
    ),
    gamlss_family = "GG",
    gamlss_optimizer = "RS",
    GAMLSS_BIC_select = FALSE,
    gamlss_continuous_X = character(0),
    gamlss_continuous_Z = character(0),
    gamlss_bic_candidates = c("linear", "pb_mu", "pb_mu_sigma"),
    gamlss_bic_trace = FALSE,
    n.cyc = 300,
    HAL_pMX = TRUE,
    HAL_pMXZ = TRUE,
    HAL_options = list(
      max_degree = 3,
      lambda_seq = exp(seq(-1, -10, length = 100)),
      num_knots = c(1000, 500, 250)
    ),
    seed = 1,
    ...
){
  if (!identical(gamlss_family, "GG")) {
    stop("The current GAMLSS density implementation supports gamlss_family = \"GG\" only.")
  }
  if (isTRUE(GAMLSS_pMX) && isTRUE(HAL_pMX)) {
    stop("Choose only one pMX density method: GAMLSS or HAL.")
  }
  if (isTRUE(GAMLSS_pMXZ) && isTRUE(HAL_pMXZ)) {
    stop("Choose only one pMXZ density method: GAMLSS or HAL.")
  }
  # number of outcomes
  p <- ncol(train_dat$Y)
  
  gamlss_selection <- list()
  
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
  
  # ----------------------------------------------------------------------
  # Build explicitly named fold-specific data frames for GAMLSS.
  # This prevents spline terms such as bs(age, ...) from being evaluated
  # against full-sample vectors in the formula environment.
  # ----------------------------------------------------------------------
  A_train <- Reduce(c, train_dat$A)
  A_valid <- Reduce(c, valid_dat$A)
  
  X_train <- as.data.frame(train_dat$X, check.names = FALSE)
  X_valid <- as.data.frame(valid_dat$X, check.names = FALSE)
  Z_train <- as.data.frame(train_dat$Z, check.names = FALSE)
  Z_valid <- as.data.frame(valid_dat$Z, check.names = FALSE)
  
  pMX_train_data <- data.frame(
    M = M_train,
    A = A_train,
    X_train,
    check.names = FALSE
  )
  
  pMX_valid_data <- data.frame(
    M = M_valid,
    A = A_valid,
    X_valid,
    check.names = FALSE
  )
  
  pMXZ_train_data <- data.frame(
    M = M_train,
    A = A_train,
    X_train,
    Z_train,
    check.names = FALSE
  )
  
  pMXZ_valid_data <- data.frame(
    M = M_valid,
    A = A_valid,
    X_valid,
    Z_valid,
    check.names = FALSE
  )
  
  if (
    nrow(pMX_train_data) != n_train ||
    nrow(pMX_valid_data) != n_valid ||
    nrow(pMXZ_train_data) != n_train ||
    nrow(pMXZ_valid_data) != n_valid
  ) {
    stop("Cross-fit GAMLSS data frames are not aligned with the fold sizes.")
  }
  
  predict_gamlss_parameters <- function(fit, newdata, fit_data) {
    list(
      mu = as.numeric(stats::predict(
        fit,
        newdata = newdata,
        data = fit_data,
        what = "mu",
        type = "response"
      )),
      sigma = as.numeric(stats::predict(
        fit,
        newdata = newdata,
        data = fit_data,
        what = "sigma",
        type = "response"
      )),
      nu = as.numeric(stats::predict(
        fit,
        newdata = newdata,
        data = fit_data,
        what = "nu",
        type = "response"
      ))
    )
  }
  
  # non-missing motion index
  index_nomissingM_train <- which(!is.na(M_train))
  index_nomissingM_valid <- which(!is.na(M_valid))
  
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
  # exposure or group indicator, conditional on baseline covariates
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
    gAn_1 <- stats::predict(gA_fit, type = "response", newdata = valid_dat$X)[[1]]
  }
  
  # fit missingness indicator regression
  if(sum(Delta_Y_train) == n_train){
    gDYn_1_AXZ_train <- gDYn_1_A0XZ_train <- gDYn_1_A1XZ_train <- rep(1, n_train)
    gDYn_1_AXZ_valid <- gDYn_1_A0XZ_valid <- gDYn_1_A1XZ_valid <- rep(1, n_valid)
  }else{
    # probability of non-missing conditioning on exposure or group indicator A, baseline covariates X and post-exposure covariates Z
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
      gDYn_1_AXZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z))[[1]]
      gDYn_1_A0XZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, train_dat$X, train_dat$Z))[[1]]
      gDYn_1_A1XZ_train <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, train_dat$X, train_dat$Z))[[1]]
      gDYn_1_AXZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z))[[1]]
      gDYn_1_A0XZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z))[[1]]
      gDYn_1_A1XZ_valid <- stats::predict(gDY_AXZ_fit, type = "response", newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z))[[1]]
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
      gDMn_1_A0_train <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, train_dat$X))[[1]]
      gDMn_1_A0_valid <- stats::predict(gDM_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X))[[1]]
    }
    
    # fit missingness indicator regression
    if(sum(Delta_Y_train) == n_train){
      gDYn_1_AX_train <- rep(1, n_train)
      gDYn_1_AX_valid <- rep(1, n_valid)
    }else{
      # probability of non-missing conditioning on exposure or group indicator A and baseline covariates X
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
        gDYn_1_AX_train <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X))[[1]]
        gDYn_1_AX_valid <- stats::predict(gDY_AX_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X))[[1]]
      }
    }
    
    ############################################
    # Motion density
    if(HAL_pMX){
      # density estimation for mediator, conditioning on exposure or group indicator and baseline covariates p(m|a,x)
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
      
      # density estimation for mediator, conditioning on exposure or group indicator and baseline covariates p(m|0,x,Delta=1)
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
    }else if(GAMLSS_pMX){
      ctrl <- gamlss::gamlss.control(n.cyc = n.cyc, c.crit = 1e-6, trace = isTRUE(gamlss_bic_trace))
      
      # GAMLSS p(M | A, X): fit on all positive observed M in the
      # training fold; evaluate the density where Delta_Y == 1.
      idx_train_all_M <- which(!is.na(M_train) & is.finite(M_train) & M_train > 0)
      idx_train_Y <- which(Delta_Y_train == 1 & !is.na(M_train) & is.finite(M_train) & M_train > 0)
      idx_valid_Y <- which(Delta_Y_valid == 1 & !is.na(M_valid) & is.finite(M_valid) & M_valid > 0)

      pMX_fit_data <- pMX_train_data[idx_train_all_M, , drop = FALSE]

      if (nrow(pMX_fit_data) == 0) {
        stop("No training observations remain for the cross-fit pMX GAMLSS model.")
      }

      pMX_spec <- .moco_make_default_gamlss_formulas(
        x_names = names(X_train),
        continuous_X = gamlss_continuous_X
      )
      pMX_fit <- .moco_fit_gamlss(
        formula = pMX_spec$mu,
        sigma.formula = pMX_spec$sigma,
        nu.formula = pMX_spec$nu,
        data = pMX_fit_data,
        control = ctrl,
        gamlss_optimizer = gamlss_optimizer,
        model_label = "cross-fit pMX GAMLSS"
      )
      gamlss_selection$pMX <- list(
        method = "default_pb_mu_sigma_nu1",
        optimizer = .moco_gamlss_optimizer_label(gamlss_optimizer),
        formulas = pMX_spec
      )

      pMXn_A_train <- rep(NA_real_, n_train)
      pMXn_A_valid <- rep(NA_real_, n_valid)
      
      if (length(idx_train_Y) > 0) {
        pMX_predict_train <- predict_gamlss_parameters(
          fit = pMX_fit,
          newdata = pMX_train_data[idx_train_Y, , drop = FALSE],
          fit_data = pMX_fit_data
        )
        
        pMXn_A_train[idx_train_Y] <- gamlss.dist::dGG(
          M_train[idx_train_Y],
          mu = pMX_predict_train$mu,
          sigma = pMX_predict_train$sigma,
          nu = pMX_predict_train$nu
        )
      }
      
      if (length(idx_valid_Y) > 0) {
        pMX_predict_valid <- predict_gamlss_parameters(
          fit = pMX_fit,
          newdata = pMX_valid_data[idx_valid_Y, , drop = FALSE],
          fit_data = pMX_fit_data
        )
        
        pMXn_A_valid[idx_valid_Y] <- gamlss.dist::dGG(
          M_valid[idx_valid_Y],
          mu = pMX_predict_valid$mu,
          sigma = pMX_predict_valid$sigma,
          nu = pMX_predict_valid$nu
        )
      }
      
      # GAMLSS p(M | A = 0, X, Delta_M = 1)
      idx_train_M <- which(Delta_M_train == 1 & !is.na(M_train) & is.finite(M_train) & M_train > 0)
      idx_valid_M <- which(Delta_M_valid == 1 & !is.na(M_valid) & is.finite(M_valid) & M_valid > 0)

      pMXD_fit_data <- pMX_train_data[idx_train_M, , drop = FALSE]

      if (nrow(pMXD_fit_data) == 0) {
        stop("No training observations remain for the cross-fit pMXD GAMLSS model.")
      }

      pMXD_spec <- .moco_make_default_gamlss_formulas(
        x_names = names(X_train),
        continuous_X = gamlss_continuous_X
      )
      pMXD_fit <- .moco_fit_gamlss(
        formula = pMXD_spec$mu,
        sigma.formula = pMXD_spec$sigma,
        nu.formula = pMXD_spec$nu,
        data = pMXD_fit_data,
        control = ctrl,
        gamlss_optimizer = gamlss_optimizer,
        model_label = "cross-fit pMXD GAMLSS"
      )
      gamlss_selection$pMXD <- list(
        method = "default_pb_mu_sigma_nu1",
        optimizer = .moco_gamlss_optimizer_label(gamlss_optimizer),
        formulas = pMXD_spec
      )

      pMXDn_A0_train <- rep(NA_real_, n_train)
      pMXDn_A0_valid <- rep(NA_real_, n_valid)
      pMXDn_A0_train[Delta_M_train == 0] <- 0
      pMXDn_A0_valid[Delta_M_valid == 0] <- 0
      
      if (length(idx_train_M) > 0) {
        pMXD_new_train <- pMX_train_data[idx_train_M, , drop = FALSE]
        pMXD_new_train$A <- 0
        
        pMXD_predict_train <- predict_gamlss_parameters(
          fit = pMXD_fit,
          newdata = pMXD_new_train,
          fit_data = pMXD_fit_data
        )
        
        pMXDn_A0_train[idx_train_M] <- gamlss.dist::dGG(
          M_train[idx_train_M],
          mu = pMXD_predict_train$mu,
          sigma = pMXD_predict_train$sigma,
          nu = pMXD_predict_train$nu
        )
      }
      
      if (length(idx_valid_M) > 0) {
        pMXD_new_valid <- pMX_valid_data[idx_valid_M, , drop = FALSE]
        pMXD_new_valid$A <- 0
        
        pMXD_predict_valid <- predict_gamlss_parameters(
          fit = pMXD_fit,
          newdata = pMXD_new_valid,
          fit_data = pMXD_fit_data
        )
        
        pMXDn_A0_valid[idx_valid_M] <- gamlss.dist::dGG(
          M_valid[idx_valid_M],
          mu = pMXD_predict_valid$mu,
          sigma = pMXD_predict_valid$sigma,
          nu = pMXD_predict_valid$nu
        )
      }
    }else{
      pMX_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X)[Delta_Y_train == 1,])
      pMXn_A_train <- rep(NA, n_train)
      pMXn_A_valid <- rep(NA, n_valid)
      pMXn_A_train[Delta_Y_train == 1] <- (1/M_train)[Delta_Y_train == 1] * dnorm(log(M_train)[Delta_Y_train == 1], mean = stats::predict(pMX_fit, newdata = data.frame(train_dat$A, train_dat$X)[Delta_Y_train == 1,]), sd = sd(pMX_fit$residuals))
      pMXn_A_valid[Delta_Y_valid == 1] <- (1/M_valid)[Delta_Y_valid == 1] * dnorm(log(M_valid)[Delta_Y_valid == 1], mean = stats::predict(pMX_fit, newdata = data.frame(valid_dat$A, valid_dat$X)[Delta_Y_valid == 1,]), sd = sd(pMX_fit$residuals))
      
      pMXD_fit <- stats::glm(paste0("log_M ~ ", glm_formula$pMX), family = gaussian(), data = data.frame(log_M = log(M_train), train_dat$A, train_dat$X)[Delta_M_train == 1,])
      pMXDn_A0_train <- rep(NA, n_train)
      pMXDn_A0_train[Delta_M_train==0] <- 0
      pMXDn_A0_train[Delta_M_train==1] <- (1/M_train)[Delta_M_train == 1] * dnorm(log(M_train)[Delta_M_train == 1], mean = stats::predict(pMXD_fit, newdata = data.frame(A=0, train_dat$X)[Delta_M_train == 1,]), sd = sd(pMXD_fit$residuals))
      pMXDn_A0_valid <- rep(NA, n_valid)
      pMXDn_A0_valid[Delta_M_valid==0] <- 0
      pMXDn_A0_valid[Delta_M_valid==1] <- (1/M_valid)[Delta_M_valid == 1] * dnorm(log(M_valid)[Delta_M_valid == 1], mean = stats::predict(pMXD_fit, newdata = data.frame(A=0, valid_dat$X)[Delta_M_valid == 1,]), sd = sd(pMXD_fit$residuals))
    }
    
    if(HAL_pMXZ){
      # density estimation for mediator, conditioning on exposure or group indicator, binary mediator,
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
      
      # density estimation for mediator, conditioning on exposure or group indicator, binary mediator,
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
    }else if(GAMLSS_pMXZ){
      ctrl <- gamlss::gamlss.control(n.cyc = n.cyc, c.crit = 1e-6, trace = isTRUE(gamlss_bic_trace))
      
      # GAMLSS p(M | A, X, Z): fit on all positive observed M in the
      # training fold; evaluate the density where Delta_Y == 1.
      idx_train_all_M <- which(!is.na(M_train) & is.finite(M_train) & M_train > 0)
      idx_train_Y <- which(Delta_Y_train == 1 & !is.na(M_train) & is.finite(M_train) & M_train > 0)
      idx_valid_Y <- which(Delta_Y_valid == 1 & !is.na(M_valid) & is.finite(M_valid) & M_valid > 0)

      pMXZ_fit_data <- pMXZ_train_data[idx_train_all_M, , drop = FALSE]

      if (nrow(pMXZ_fit_data) == 0) {
        stop("No training observations remain for the cross-fit pMXZ GAMLSS model.")
      }

      pMXZ_spec <- .moco_make_default_gamlss_formulas(
        x_names = names(X_train),
        z_names = names(Z_train),
        continuous_X = gamlss_continuous_X,
        continuous_Z = gamlss_continuous_Z
      )
      pMXZ_fit <- .moco_fit_gamlss(
        formula = pMXZ_spec$mu,
        sigma.formula = pMXZ_spec$sigma,
        nu.formula = pMXZ_spec$nu,
        data = pMXZ_fit_data,
        control = ctrl,
        gamlss_optimizer = gamlss_optimizer,
        model_label = "cross-fit pMXZ GAMLSS"
      )
      gamlss_selection$pMXZ <- list(
        method = "default_pb_mu_sigma_nu1",
        optimizer = .moco_gamlss_optimizer_label(gamlss_optimizer),
        formulas = pMXZ_spec
      )

      pMXZn_A_train <- rep(NA_real_, n_train)
      pMXZn_A0_train <- rep(NA_real_, n_train)
      pMXZn_A1_train <- rep(NA_real_, n_train)
      
      pMXZn_A_valid <- rep(NA_real_, n_valid)
      pMXZn_A0_valid <- rep(NA_real_, n_valid)
      pMXZn_A1_valid <- rep(NA_real_, n_valid)
      
      if (length(idx_train_Y) > 0) {
        train_obs <- pMXZ_train_data[idx_train_Y, , drop = FALSE]
        train_A0 <- train_obs
        train_A1 <- train_obs
        train_A0$A <- 0
        train_A1$A <- 1
        
        pMXZ_predict_train <- predict_gamlss_parameters(
          fit = pMXZ_fit,
          newdata = train_obs,
          fit_data = pMXZ_fit_data
        )
        pMXZ_predict0_train <- predict_gamlss_parameters(
          fit = pMXZ_fit,
          newdata = train_A0,
          fit_data = pMXZ_fit_data
        )
        pMXZ_predict1_train <- predict_gamlss_parameters(
          fit = pMXZ_fit,
          newdata = train_A1,
          fit_data = pMXZ_fit_data
        )
        
        pMXZn_A_train[idx_train_Y] <- gamlss.dist::dGG(
          M_train[idx_train_Y],
          mu = pMXZ_predict_train$mu,
          sigma = pMXZ_predict_train$sigma,
          nu = pMXZ_predict_train$nu
        )
        pMXZn_A0_train[idx_train_Y] <- gamlss.dist::dGG(
          M_train[idx_train_Y],
          mu = pMXZ_predict0_train$mu,
          sigma = pMXZ_predict0_train$sigma,
          nu = pMXZ_predict0_train$nu
        )
        pMXZn_A1_train[idx_train_Y] <- gamlss.dist::dGG(
          M_train[idx_train_Y],
          mu = pMXZ_predict1_train$mu,
          sigma = pMXZ_predict1_train$sigma,
          nu = pMXZ_predict1_train$nu
        )
      }
      
      if (length(idx_valid_Y) > 0) {
        valid_obs <- pMXZ_valid_data[idx_valid_Y, , drop = FALSE]
        valid_A0 <- valid_obs
        valid_A1 <- valid_obs
        valid_A0$A <- 0
        valid_A1$A <- 1
        
        pMXZ_predict_valid <- predict_gamlss_parameters(
          fit = pMXZ_fit,
          newdata = valid_obs,
          fit_data = pMXZ_fit_data
        )
        pMXZ_predict0_valid <- predict_gamlss_parameters(
          fit = pMXZ_fit,
          newdata = valid_A0,
          fit_data = pMXZ_fit_data
        )
        pMXZ_predict1_valid <- predict_gamlss_parameters(
          fit = pMXZ_fit,
          newdata = valid_A1,
          fit_data = pMXZ_fit_data
        )
        
        pMXZn_A_valid[idx_valid_Y] <- gamlss.dist::dGG(
          M_valid[idx_valid_Y],
          mu = pMXZ_predict_valid$mu,
          sigma = pMXZ_predict_valid$sigma,
          nu = pMXZ_predict_valid$nu
        )
        pMXZn_A0_valid[idx_valid_Y] <- gamlss.dist::dGG(
          M_valid[idx_valid_Y],
          mu = pMXZ_predict0_valid$mu,
          sigma = pMXZ_predict0_valid$sigma,
          nu = pMXZ_predict0_valid$nu
        )
        pMXZn_A1_valid[idx_valid_Y] <- gamlss.dist::dGG(
          M_valid[idx_valid_Y],
          mu = pMXZ_predict1_valid$mu,
          sigma = pMXZ_predict1_valid$sigma,
          nu = pMXZ_predict1_valid$nu
        )
      }
      
      # GAMLSS p(M | A, X, Z, Delta_M = 1)
      idx_train_M <- which(Delta_M_train == 1 & !is.na(M_train) & is.finite(M_train) & M_train > 0)
      idx_valid_M <- which(Delta_M_valid == 1 & !is.na(M_valid) & is.finite(M_valid) & M_valid > 0)

      pMXZD_fit_data <- pMXZ_train_data[idx_train_M, , drop = FALSE]

      if (nrow(pMXZD_fit_data) == 0) {
        stop("No training observations remain for the cross-fit pMXZD GAMLSS model.")
      }

      pMXZD_spec <- .moco_make_default_gamlss_formulas(
        x_names = names(X_train),
        z_names = names(Z_train),
        continuous_X = gamlss_continuous_X,
        continuous_Z = gamlss_continuous_Z
      )
      pMXZD_fit <- .moco_fit_gamlss(
        formula = pMXZD_spec$mu,
        sigma.formula = pMXZD_spec$sigma,
        nu.formula = pMXZD_spec$nu,
        data = pMXZD_fit_data,
        control = ctrl,
        gamlss_optimizer = gamlss_optimizer,
        model_label = "cross-fit pMXZD GAMLSS"
      )
      gamlss_selection$pMXZD <- list(
        method = "default_pb_mu_sigma_nu1",
        optimizer = .moco_gamlss_optimizer_label(gamlss_optimizer),
        formulas = pMXZD_spec
      )

      pMXZDn_A_train <- rep(NA_real_, n_train)
      pMXZDn_A_valid <- rep(NA_real_, n_valid)
      pMXZDn_A_train[Delta_M_train == 0] <- 0
      pMXZDn_A_valid[Delta_M_valid == 0] <- 0
      
      if (length(idx_train_M) > 0) {
        pMXZD_predict_train <- predict_gamlss_parameters(
          fit = pMXZD_fit,
          newdata = pMXZ_train_data[idx_train_M, , drop = FALSE],
          fit_data = pMXZD_fit_data
        )
        
        pMXZDn_A_train[idx_train_M] <- gamlss.dist::dGG(
          M_train[idx_train_M],
          mu = pMXZD_predict_train$mu,
          sigma = pMXZD_predict_train$sigma,
          nu = pMXZD_predict_train$nu
        )
      }
      
      if (length(idx_valid_M) > 0) {
        pMXZD_predict_valid <- predict_gamlss_parameters(
          fit = pMXZD_fit,
          newdata = pMXZ_valid_data[idx_valid_M, , drop = FALSE],
          fit_data = pMXZD_fit_data
        )
        
        pMXZDn_A_valid[idx_valid_M] <- gamlss.dist::dGG(
          M_valid[idx_valid_M],
          mu = pMXZD_predict_valid$mu,
          sigma = pMXZD_predict_valid$sigma,
          nu = pMXZD_predict_valid$nu
        )
      }
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
      pMXZDn_A_train <- rep(NA, n_train)
      pMXZDn_A_train[Delta_M_train==0] <- 0
      pMXZDn_A_train[Delta_M_train==1] <- (1/M_train)[Delta_M_train==1] * dnorm(log(M_train)[Delta_M_train==1], mean = stats::predict(pMXZD_fit, newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z)[Delta_M_train==1,]),
                                                                                sd = sd(pMXZD_fit$residuals))
      pMXZDn_A_valid <- rep(NA, n_valid)
      pMXZDn_A_valid[Delta_M_valid==0] <- 0
      pMXZDn_A_valid[Delta_M_valid==1] <- (1/M_valid)[Delta_M_valid==1] * dnorm(log(M_valid)[Delta_M_valid==1], mean = stats::predict(pMXZD_fit, newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z)[Delta_M_valid==1,]),
                                                                                sd = sd(pMXZD_fit$residuals))
    }
  }
  
  # define parameters for storing results
  est <- matrix(nrow = 2, ncol = p)
  eif_mat <- list(p)
  cov_mat <- list(p)
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
      mu_MXZn_A0_train[index_nomissingM_train] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, train_dat$M, train_dat$X, train_dat$Z)[index_nomissingM_train,])[[1]]
      mu_MXZn_A1_train[index_nomissingM_train] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, train_dat$M, train_dat$X, train_dat$Z)[index_nomissingM_train,])[[1]]
      mu_MXZn_A0_valid[index_nomissingM_valid] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 0, valid_dat$M, valid_dat$X, valid_dat$Z)[index_nomissingM_valid,])[[1]]
      mu_MXZn_A1_valid[index_nomissingM_valid] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(A = 1, valid_dat$M, valid_dat$X, valid_dat$Z)[index_nomissingM_valid,])[[1]]
      mu_MXZn_A_train[index_nomissingM_train] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(train_dat$A, train_dat$M, train_dat$X, train_dat$Z)[index_nomissingM_train,])[[1]]
      mu_MXZn_A_valid[index_nomissingM_valid] <- stats::predict(mu_AMXZ_fit, newdata = data.frame(valid_dat$A, valid_dat$M, valid_dat$X, valid_dat$Z)[index_nomissingM_valid,])[[1]]
    }
    
    if(M_indicator){
      # fit additional pseudo-outcome regression for QY*qM/qMZ, conditioning on exposure or group indicator, binary mediator and baseline covariates
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
        eta_AXZn_A_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(train_dat$A, train_dat$X, train_dat$Z))[[1]]
        eta_AXZn_A0_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, train_dat$X, train_dat$Z))[[1]]
        eta_AXZn_A1_train <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, train_dat$X, train_dat$Z))[[1]]
        
        eta_AXZn_A_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(valid_dat$A, valid_dat$X, valid_dat$Z))[[1]]
        eta_AXZn_A0_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 0, valid_dat$X, valid_dat$Z))[[1]]
        eta_AXZn_A1_valid <- stats::predict(eta_AXZ_fit, type = "response", newdata = data.frame(A = 1, valid_dat$X, valid_dat$Z))[[1]]
      }
      
      # fit pseudo-outcome regression for Qd, conditioning on exposure or group indicator, mediator and baseline covariates
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
        eta_AXMn_A0_train[index_nomissingM_train] <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, train_dat$M, train_dat$X)[index_nomissingM_train,])[[1]]
        eta_AXMn_A1_train[index_nomissingM_train] <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, train_dat$M, train_dat$X)[index_nomissingM_train,])[[1]]
        
        eta_AXMn_A0_valid[index_nomissingM_valid] <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 0, valid_dat$M, valid_dat$X)[index_nomissingM_valid,])[[1]]
        eta_AXMn_A1_valid[index_nomissingM_valid] <- stats::predict(eta_AXM_fit, newdata = data.frame(A = 1, valid_dat$M, valid_dat$X)[index_nomissingM_valid,])[[1]]
      }else{
        set.seed(seed)
        eta_AXM_fit <- hal9001::fit_hal(X = data.frame(train_dat$A, train_dat$M, train_dat$X)[Delta_Y_train==1,], Y = mu_pseudo_A_star_train[Delta_Y_train==1])
        eta_AXMn_A0_train <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 0, train_dat$M, train_dat$X))
        eta_AXMn_A1_train <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 1, train_dat$M, train_dat$X))
        
        eta_AXMn_A0_valid <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 0, valid_dat$M, valid_dat$X))
        eta_AXMn_A1_valid <- stats::predict(eta_AXM_fit, new_data = data.frame(A = 1, valid_dat$M, valid_dat$X))
      }
      
    }else{
      eta_AXZn_A_train = mu_MXZn_A_train
    }
    
    # fit pseudo-outcome regression for u_star, conditioning on exposure or group indicator and baseline covariates
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
      xi_AXn_A0 <- stats::predict(xi_fit, newdata = data.frame(A = 0, valid_dat$X))[[1]]
      xi_AXn_A1 <- stats::predict(xi_fit, newdata = data.frame(A = 1, valid_dat$X))[[1]]
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
  out <- list(
    est = est,
    eif_mat = eif_mat,
    gamlss_selection = gamlss_selection
  )
  
  return(out)
}

one_step_cross_gamlss <- function(
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
    GAMLSS_pMX = FALSE,
    GAMLSS_pMXZ = FALSE,
    gamlss_formula = list(
      pMX_mu = NULL,
      pMX_sigma = ~ 1,
      pMX_nu = ~ 1,
      pMXZ_mu = NULL,
      pMXZ_sigma = ~ 1,
      pMXZ_nu = ~ 1
    ),
    gamlss_family = "GG",
    gamlss_optimizer = "RS",
    GAMLSS_BIC_select = FALSE,
    gamlss_continuous_X = character(0),
    gamlss_continuous_Z = character(0),
    gamlss_bic_candidates = c("linear", "pb_mu", "pb_mu_sigma"),
    gamlss_bic_trace = FALSE,
    n.cyc = 300,
    HAL_pMX = TRUE,
    HAL_pMXZ = TRUE,
    HAL_options = list(max_degree = 3, lambda_seq = exp(seq(-1, -10, length = 100)), num_knots = c(1000, 500, 250)),
    seed = 1,
    cv_folds = 5,
    ...
){
  
  #from moco
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
  
  # number of participants
  n <- nrow(Y)
  # number of outcomes
  p <- ncol(Y)
  # whether motion is available
  M_indicator <- (!is.null(M))
  
  if(!M_indicator){
    # missingness indicator as in Nebel
    Delta_Y <- ifelse(Delta_M + Delta_Y == 2, 1, 0)
    # impute M to be all 1's for computation convenience
    M <- rep(1, n)
    # M[Delta_Y == 0] <- NA
  }else if(!is.null(thresh)){
    # if thresh is not null, collapse M into dummy variables based on truncated level thresh
    Delta_M <- as.numeric(M < thresh)
    Delta_M[is.na(M)] <- 0
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
  eif_mat <- list(p)
  gamlss_selection_by_fold <- vector("list", cv_folds)
  for(j in 1:p){
    eif_mat[[j]] <- matrix(nrow = n, ncol = 2)
  }
  
  # cross fitting
  for(i in 1:cv_folds){
    rslt_cross <- fit_mechanism_gamlss(
      train_dat = lapply(dat, function(x){x[-folds[[i]], , drop = FALSE]}),
      valid_dat = lapply(dat, function(x){x[folds[[i]], , drop = FALSE]}),
      M_indicator = M_indicator,
      SL_library = SL_library,
      SL_library_customize = SL_library_customize,
      glm_formula = glm_formula,
      GAMLSS_pMX = GAMLSS_pMX,
      GAMLSS_pMXZ = GAMLSS_pMXZ,
      gamlss_formula = gamlss_formula,
      gamlss_family = gamlss_family,
      gamlss_optimizer = gamlss_optimizer,
      GAMLSS_BIC_select = GAMLSS_BIC_select,
      gamlss_continuous_X = gamlss_continuous_X,
      gamlss_continuous_Z = gamlss_continuous_Z,
      gamlss_bic_candidates = gamlss_bic_candidates,
      gamlss_bic_trace = gamlss_bic_trace,
      n.cyc = n.cyc,
      HAL_pMX = HAL_pMX,
      HAL_pMXZ = HAL_pMXZ,
      HAL_options = HAL_options,
      seed = seed
    )
    
    gamlss_selection_by_fold[[i]] <- rslt_cross$gamlss_selection
    
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
  out <- list(
    est = est,
    adj_association = adj_association,
    eif_mat = eif_mat,
    cov_mat = cov_mat,
    gamlss_selection_by_fold = gamlss_selection_by_fold
  )
}

