# GAMLSS (generalized gamma) conditional densities of motion for MoCo.
#
# These helpers are used when density_method = "gamlss". They fit
#   pMX, pMXZ   : p(m | a, x[, z], Delta_Y = 1), GG on the Delta_Y = 1 participants
#   pMXD, pMXZD : p(m | a, x[, z], Delta_M = 1), GG on the Delta_M = 1 participants.
#

# ---------------------------------------------------------------------------
# options
# ---------------------------------------------------------------------------

# fixed fitting settings of gamlss::gamlss() (not user options)
gg_fit_control <- list(method = "RS", n.cyc = 500, c.crit = 0.001)

gg_default_options <- function(){
  list(
    mu = "smooth",
    sigma = "linear",
    nu = "1",
    smooth = NULL,
    trace = TRUE
  )
}

# fill in defaults and check the user's gamlss_options
gg_check_options <- function(gamlss_options){
  opts <- gg_default_options()
  
  if(!is.null(gamlss_options)){
    if(!is.list(gamlss_options)){
      stop("gamlss_options must be a list.")
    }
    
    unknown <- setdiff(names(gamlss_options), names(opts))
    if(length(unknown) > 0){
      stop("Unknown entries in gamlss_options: ", paste(unknown, collapse = ", "),
           ". Allowed entries: ", paste(names(opts), collapse = ", "), ".")
    }
    
    for(nm in names(gamlss_options)){
      if(!is.null(gamlss_options[[nm]]) || nm == "smooth"){
        opts[nm] <- list(gamlss_options[[nm]])
      }
    }
  }
  
  if(!(is.character(opts$mu) && length(opts$mu) == 1 && opts$mu %in% c("linear", "smooth")) &&
     !inherits(opts$mu, "formula")){
    stop("gamlss_options$mu must be \"linear\", \"smooth\", or a one-sided formula.")
  }
  if(!(is.character(opts$sigma) && length(opts$sigma) == 1 && opts$sigma %in% c("linear", "smooth", "1", "A")) &&
     !inherits(opts$sigma, "formula")){
    stop("gamlss_options$sigma must be \"linear\", \"smooth\", \"1\", \"A\", or a one-sided formula.")
  }
  if(identical(opts$nu, 1)) opts$nu <- "1"
  if(identical(opts$nu, "0")) opts$nu <- 0
  if(!(identical(opts$nu, "1") || identical(as.numeric(opts$nu), 0))){
    stop("gamlss_options$nu must be \"1\" (nu estimated, constant) or 0 (nu fixed at 0, the log-normal).")
  }
  if(!is.null(opts$smooth) && !is.character(opts$smooth)){
    stop("gamlss_options$smooth must be NULL or a character vector of continuous covariate names.")
  }
  if(!(is.logical(opts$trace) && length(opts$trace) == 1 && !is.na(opts$trace))){
    stop("gamlss_options$trace must be TRUE or FALSE.")
  }
  
  opts
}

# check that the suggested packages are available
gg_check_packages <- function(){
  for(pkg in c("gamlss", "gamlss.dist")){
    if(!requireNamespace(pkg, quietly = TRUE)){
      stop("Package '", pkg, "' is required for density_method = \"gamlss\". ",
           "Please install it with install.packages(\"", pkg, "\").")
    }
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# formulas
# ---------------------------------------------------------------------------

# plain data frame: drop AsIs class, keep syntactic names
gg_clean_data <- function(W){
  W <- as.data.frame(W)
  for(j in seq_along(W)){
    if(inherits(W[[j]], "AsIs")){
      W[[j]] <- unclass(W[[j]])
    }
  }
  W
}

# covariates that enter the model: non-constant columns
gg_covariates <- function(W){
  keep <- vapply(W, function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
  names(W)[keep]
}

# continuous covariates that receive pb(): the user's list, or numeric columns with > 10 distinct values
gg_continuous <- function(W, covars, smooth){
  if(!is.null(smooth)){
    return(intersect(smooth, setdiff(covars, "A")))
  }
  cont <- vapply(covars, function(v){
    v != "A" && is.numeric(W[[v]]) && length(unique(W[[v]])) > 10
  }, logical(1))
  covars[cont]
}

gg_rhs <- function(type, covars, cont){
  if(identical(type, "1") || length(covars) == 0){
    return("1")
  }
  if(identical(type, "A")){
    return(if("A" %in% covars) "A" else "1")
  }
  terms <- covars
  if(identical(type, "smooth")){
    terms <- ifelse(covars %in% cont, paste0("pb(", covars, ")"), covars)
  }
  paste(terms, collapse = " + ")
}

# mu, sigma and nu formulas for one density component
gg_formulas <- function(W, label, opts){
  covars <- gg_covariates(W)
  cont <- gg_continuous(W, covars, opts$smooth)

  as_rhs <- function(f){
    if(inherits(f, "formula")) paste(deparse(f[[length(f)]], width.cutoff = 500L), collapse = " ") else f
  }
  mu_rhs <- if(inherits(opts$mu, "formula")) as_rhs(opts$mu) else gg_rhs(opts$mu, covars, cont)
  sigma_rhs <- if(inherits(opts$sigma, "formula")) as_rhs(opts$sigma) else gg_rhs(opts$sigma, covars, cont)
  nu_fixed <- identical(suppressWarnings(as.numeric(opts$nu)), 0)

  # every variable used must be in W
  used <- unique(c(all.vars(stats::as.formula(paste("~", mu_rhs))),
                   all.vars(stats::as.formula(paste("~", sigma_rhs)))))
  missing_vars <- setdiff(used, names(W))
  if(length(missing_vars) > 0){
    stop("gamlss formula for ", label, " uses variables that are not available for this density: ",
         paste(missing_vars, collapse = ", "), ". Note that pMX and pMXD only use A and X.")
  }
  list(mu = paste(".motion ~", mu_rhs), 
       sigma = paste("~", sigma_rhs), 
       nu = "~ 1",
       nu_fixed = nu_fixed)
}

# ---------------------------------------------------------------------------
# fitting
# ---------------------------------------------------------------------------

# fit one density component with gamlss::gamlss() and the GG family.
# There is no fallback: if the fit fails or does not converge, moco() stops with an error.
gg_component <- function(m, W, label, opts){
  W <- gg_clean_data(W)
  dat <- data.frame(.motion = m, W)
  forms <- gg_formulas(W, label, opts)
  env <- new.env(parent = asNamespace("gamlss"))
  env$.dat <- dat
  mk <- function(txt){ f <- stats::as.formula(txt); environment(f) <- env; f }
  args <- list(formula = mk(forms$mu), 
               sigma.formula = mk(forms$sigma), 
               nu.formula = mk(forms$nu),
               family = quote(gamlss.dist::GG()), 
               data = quote(.dat), 
               method = call(gg_fit_control$method),
               control = gamlss::gamlss.control(n.cyc = gg_fit_control$n.cyc, 
                                                c.crit = gg_fit_control$c.crit,
                                                trace = opts$trace))
  if(forms$nu_fixed){
    args$nu.fix <- TRUE
    args$nu.start <- 0
  }
  
  cl <- as.call(c(list(as.name("gamlss")), args))
  fit <- tryCatch(suppressWarnings(eval(cl, env)), error = function(e) conditionMessage(e))
  
  if(!inherits(fit, "gamlss") || !isTRUE(fit$converged) || !is.finite(fit$G.deviance)){
    stop("The GAMLSS fit of ", label, " (n = ", nrow(dat), ") did not converge",
         if(is.character(fit)) paste0(": ", substr(fit, 1, 200)) else ".",
         " Try a simpler model, e.g. gamlss_options = list(mu = \"linear\", sigma = \"A\"), ",
         "or density_method = \"glm\" or \"hal\".", call. = FALSE)
  }
  list(fit = fit, 
       data = dat, 
       range_M = range(m), 
       nu_fixed = forms$nu_fixed)
}

# evaluate a fitted component at motion values m and covariates W.
# zero_outside = TRUE sets the density to 0 outside the training range of M
# (used for pMXD/pMXZD, the same rule as the GLM and HAL densities).
# thresh = c (pMXD/pMXZD with Delta_M = 1(M <= c)): obj is an untruncated GG fitted on
# all participants with observed M, and the density is truncated analytically,
# p(m | Delta_M = 1, w) = dGG(m | w) / pGG(c | w) for m <= c and 0 for m > c,
# with pGG(c | w) floored at floor_cdf.
gg_density <- function(obj, m, W, zero_outside = FALSE, thresh = NULL, floor_cdf = 1e-3){
  if(length(m) == 0){
    return(numeric(0))
  }
  W <- gg_clean_data(W)
  newdat <- data.frame(.motion = m, W)
  pars <- lapply(c(mu = "mu", sigma = "sigma", nu = "nu"), function(par){
    if(par == "nu" && isTRUE(obj$nu_fixed)){
      return(rep(0, length(m)))   # nu fixed at 0 has no nu formula to predict from
    }
    suppressWarnings(as.numeric(stats::predict(obj$fit, what = par, newdata = newdat, data = obj$data, type = "response")))
  })
  d <- gamlss.dist::dGG(m, mu = pars$mu, sigma = pars$sigma, nu = pars$nu)
  if(!is.null(thresh)){
    Fc <- gamlss.dist::pGG(thresh, mu = pars$mu, sigma = pars$sigma, nu = pars$nu)
    if(mean(Fc < floor_cdf) > 0.01){
      warning("For more than 1% of participants, the fitted P(M <= thresh) is below ", floor_cdf,
              ", so the truncated GAMLSS density uses that floor. The motion density model may fit poorly near thresh.",
              call. = FALSE)
    }
    d <- d / pmax(Fc, floor_cdf)
    d[m > thresh] <- 0
    return(d)
  }
  if(zero_outside){
    d[m < obj$range_M[1] | m > obj$range_M[2]] <- 0
  }
  d
}
