# GAMLSS (generalized gamma) conditional densities of motion for MoCo.
#
# These helpers are used when density_method = "gamlss". They fit
#   pMX, pMXZ   : p(m | a, x[, z], Delta_Y = 1), GG on the Delta_Y = 1 participants
#   pMXD, pMXZD : p(m | a, x[, z], Delta_M = 1), GG on the Delta_M = 1 participants.
#

# ---------------------------------------------------------------------------
# options
# ---------------------------------------------------------------------------
#' Default GG motion-density options
#'
#' Construct the default GAMLSS option list.
#'
#' @return A list with mu, sigma, nu, smooth and trace.
#' @noRd
gg_default_options <- function(){
  list(
    mu = "smooth",
    sigma = "linear",
    nu = "1",
    smooth = NULL,
    method = "RS",
    n.cyc = 500,
    c.crit = 0.001,
    trace = TRUE
  )
}


#' Validate GG motion-density options
#' 
#' @details Combine user settings with defaults and check supported values. nu = "1" means an estimated constant shape; nu = 0 fixes the shape at the lognormal limit.
#'
#' @param gamlss_options NULL or a list of GAMLSS density options; see moco().
#' @return The completed option list.
#' @noRd
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
  if(!(is.character(opts$method) && length(opts$method) == 1 &&
       !is.na(opts$method) && opts$method %in% c("RS", "CG", "mixed"))){
    stop('gamlss_options$method must be "RS", "CG", or "mixed".')
  }
  if(!(is.numeric(opts$n.cyc) && length(opts$n.cyc) == 1 &&
       is.finite(opts$n.cyc) && opts$n.cyc >= 1 &&
       opts$n.cyc == floor(opts$n.cyc))){
    stop("gamlss_options$n.cyc must be a positive whole number.")
  }
  if(!(is.numeric(opts$c.crit) && length(opts$c.crit) == 1 &&
       is.finite(opts$c.crit) && opts$c.crit > 0)){
    stop("gamlss_options$c.crit must be a positive finite number.")
  }
  if(!(is.logical(opts$trace) && length(opts$trace) == 1 && !is.na(opts$trace))){
    stop("gamlss_options$trace must be TRUE or FALSE.")
  }
  
  return(opts)
}

#' Check GG density dependencies
#'
#' Check availability of gamlss and gamlss.dist before fitting.
#'
#' @return TRUE invisibly, or an error naming a missing package.
#' @noRd
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

#' Prepare GAMLSS covariate data
#'
#' Convert covariates to a data frame and remove column-level AsIs classes.
#'
#' @param W Covariate data frame or an object coercible to one.
#' @return A data frame with the same covariate values.
#' @noRd
gg_clean_data <- function(W){
  W <- as.data.frame(W)
  for(j in seq_along(W)){
    if(inherits(W[[j]], "AsIs")){
      W[[j]] <- unclass(W[[j]])
    }
  }
  W
}

#' Identify nonconstant covariates
#'
#' Select columns with more than one distinct nonmissing value.
#'
#' @param W Covariate data frame.
#' @return A character vector of selected column names.
#' @noRd
gg_covariates <- function(W){
  keep <- vapply(W, function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
  names(W)[keep]
}

#' Select covariates for P-spline terms
#'
#' Choose variables to wrap in pb() when a smooth model is requested.
#'
#' @details Automatic selection uses numeric variables with more than 10
#' distinct values. A is never smoothed.
#'
#' @param W Covariate data frame.
#' @param covars Names of eligible nonconstant covariates.
#' @param smooth NULL for automatic selection, or explicit covariate names.
#' @return A character vector of selected names, excluding A.
#' @noRd
gg_continuous <- function(W, covars, smooth){
  if(!is.null(smooth)){
    return(intersect(smooth, setdiff(covars, "A")))
  }
  cont <- vapply(covars, function(v){
    v != "A" && is.numeric(W[[v]]) && length(unique(W[[v]])) > 10
  }, logical(1))
  covars[cont]
}

#' Build a GAMLSS formula right-hand side
#'
#' Construct constant, group-only, linear or smooth terms.
#'
#' @param type One of "1", "A", "linear" or "smooth".
#' @param covars Names of nonconstant covariates.
#' @param cont Names to wrap in pb() for a smooth specification.
#' @return One character string containing a formula right-hand side.
#' @noRd
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

#' Build the three GG parameter specifications
#'
#' Construct location and scale formulas and identify the fixed-zero shape option.
#'
#' @param W Covariate data frame for this density component.
#' @param label Component label used in error messages.
#' @param opts Validated list returned by gg_check_options().
#' @return A list with mu, sigma and nu formula strings and a nu_fixed flag.
#' @noRd
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

#' Fit one GG motion-density component
#'
#' Fit an ordinary parent or selected-sample GG model with gamlss().
#'
#' @details The caller chooses the training rows. With a known motion
#' threshold, pMXD/pMXZD callers supply all observed-motion training rows
#' to fit a parent model. pMX/pMXZ callers use Delta_Y = 1 rows. RS uses
#' at most 500 cycles. trace follows opts$trace.
#'
#' @param m Positive observed training motion values.
#' @param W Aligned training covariates: A and X, optionally Z.
#' @param label Component name: pMX, pMXD, pMXZ or pMXZD.
#' @param opts Validated list returned by gg_check_options().
#' @return A list with fit, data, range_M and nu_fixed.
#' @noRd
gg_component <- function(m, W, label, opts){
  W <- gg_clean_data(W)
  dat <- data.frame(.motion = m, W)
  forms <- gg_formulas(W, label, opts)
  env <- new.env(parent = asNamespace("gamlss"))
  env$.dat <- dat
  mk <- function(txt){ f <- stats::as.formula(txt); environment(f) <- env; f }
  # These calls are evaluated inside gamlss(), where its algorithms are defined.
  # mixed() keeps GAMLSS defaults: 1 RS iteration, then up to 20 CG iterations.
  # Its iteration schedule is independent of control$n.cyc.
  method_call <- call(opts$method)
  args <- list(formula = mk(forms$mu), 
               sigma.formula = mk(forms$sigma), 
               nu.formula = mk(forms$nu),
               family = quote(gamlss.dist::GG()), 
               data = quote(.dat), 
               method = method_call,
               control = gamlss::gamlss.control(n.cyc = opts$n.cyc, 
                                                c.crit = opts$c.crit,
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
  return(list(fit = fit, 
       data = dat, 
       range_M = range(m), 
       nu_fixed = forms$nu_fixed)
  )
}

#' Evaluate the GG conditioning CDF
#'
#' Use the lognormal limit for zero shape and pGG() otherwise.
#'
#' @details For abs(nu) <= 1e-6, evaluate plnorm(q, log(mu), sigma),
#' avoiding the vector pGG() zero-shape bug. Current supported models
#' have constant nu; mixed-sign vector inputs are not a supported
#' prediction path in this implementation.
#'
#' @param q Evaluation cutoff or vector of cutoffs.
#' @param mu Positive GG location parameter(s).
#' @param sigma Positive GG scale parameter(s).
#' @param nu GG shape parameter(s); current callers supply a constant shape per model.
#' @return A numeric vector of lower-tail probabilities.
#' @noRd
gg_cdf <- function(q, mu, sigma, nu){
  n <- max(length(q), length(mu), length(sigma), length(nu))
  q <- rep_len(q, n); mu <- rep_len(mu, n); sigma <- rep_len(sigma, n); nu <- rep_len(nu, n)
  out <- numeric(n)
  ln <- abs(nu) <= 1e-6
  if(any(ln)) out[ln] <- stats::plnorm(q[ln], meanlog = log(mu[ln]), sdlog = sigma[ln])
  if(any(!ln)) out[!ln] <- gamlss.dist::pGG(q[!ln], mu = mu[!ln], sigma = sigma[!ln], nu = nu[!ln])
  
  return(out)
}

#' Evaluate a GG motion-density component
#'
#' Predict GG parameters and optionally condition the density at a motion cutoff.
#'
#' @details With thresh, obj must be a parent fit: return f(m)/F(thresh)
#' for m <= thresh and zero above it. This branch bypasses empirical
#' range clipping. Activating floor_cdf regularizes the result and changes
#' its total mass. mu and sigma predictions use their response scales.
#'
#' @param obj Fitted component returned by gg_component().
#' @param m Observed positive motion values at which to evaluate the density.
#' @param W Prediction covariates aligned to m.
#' @param zero_outside If TRUE and thresh is NULL, zero predictions outside the training motion range.
#' @param thresh NULL, or a positive cutoff for analytical conditioning on M <= thresh.
#' @param floor_cdf Lower bound for the conditioning CDF; default 0.001.
#' @return A numeric vector of density predictions.
#' @noRd
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
    Fc <- gg_cdf(thresh, mu = pars$mu, sigma = pars$sigma, nu = pars$nu)
    if(mean(Fc < floor_cdf) > 0.01){
      warning("For more than 1% of participants, the fitted P(M <= thresh) is below ", floor_cdf,
              ", so the truncated GAMLSS density uses that floor.",
              call. = FALSE)
    }
    d <- d / pmax(Fc, floor_cdf)
    d[m > thresh] <- 0
    return(d)
  }
  if(zero_outside){
    d[m < obj$range_M[1] | m > obj$range_M[2]] <- 0
  }
  
  return(d)
}
