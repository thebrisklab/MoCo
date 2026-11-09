# GAMLSS (generalized gamma) conditional densities of motion for MoCo.
#
# These helpers are used when density_method = "gamlss". They fit
#   pMX, pMXZ   : p(m | a, x[, z], Delta_Y = 1), GG on the Delta_Y = 1 participants
#   pMXD, pMXZD : p(m | a, x[, z], Delta_M = 1), GG on the Delta_M = 1 participants,
#                 right-truncated at Delta_M_threshold when Delta_M = 1(M < threshold).
#
# Naming note: gamlss::pb() finds its working frame by searching the call stack for
# a call whose text starts with "gamlss". No helper here may have a name that starts
# with "gamlss", otherwise pb() would pick the wrong frame.

# ---------------------------------------------------------------------------
# options
# ---------------------------------------------------------------------------

gg_default_options <- function(){
  list(
    family = "GG",
    mu = "smooth",
    sigma = "linear",
    nu = "1",
    smooth = NULL,
    formula = NULL,
    method = c("RS", "CG", "mixed"),
    n.cyc = 500,
    c.crit = 0.001,
    nu0_fallback = TRUE,
    trace = FALSE
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
      if(!is.null(gamlss_options[[nm]]) || nm %in% c("smooth", "formula")){
        opts[nm] <- list(gamlss_options[[nm]])
      }
    }
  }
  if(!identical(opts$family, "GG")){
    stop("gamlss_options$family: only \"GG\" (generalized gamma) is currently supported.")
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
  if(!all(opts$method %in% c("RS", "CG", "mixed")) || length(opts$method) == 0){
    stop("gamlss_options$method must contain one or more of \"RS\", \"CG\", \"mixed\".")
  }
  if(!is.null(opts$formula)){
    if(!is.list(opts$formula) || !all(names(opts$formula) %in% c("pMX", "pMXD", "pMXZ", "pMXZD"))){
      stop("gamlss_options$formula must be a list with entries among pMX, pMXD, pMXZ, pMXZD, each a list(mu =, sigma =, nu =).")
    }
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

# resolve density_method into one method for pMX/pMXD and one for pMXZ/pMXZD
resolve_density_method <- function(density_method){
  if(is.null(density_method)){
    density_method <- "hal"
  }
  nms <- names(density_method)
  if(!is.character(density_method) || !(length(density_method) %in% c(1, 2)) ||
     !all(tolower(density_method) %in% c("glm", "hal", "gamlss"))){
    stop("density_method must be \"hal\", \"glm\" or \"gamlss\", or a length-2 vector c(pMX = , pMXZ = ) of these.")
  }
  density_method <- tolower(density_method)
  names(density_method) <- nms
  if(length(density_method) == 1){
    density_method <- rep(density_method, 2)
  }
  if(!is.null(names(density_method)) && all(names(density_method) %in% c("pMX", "pMXZ"))){
    density_method <- density_method[c("pMX", "pMXZ")]
  }
  names(density_method) <- c("pMX", "pMXZ")
  density_method
}

# fill in and check the right-hand sides of the log-normal GLM motion densities
check_glm_options <- function(glm_options){
  out <- list(pMX = ".", pMXZ = ".")
  if(is.null(glm_options)){
    return(out)
  }
  if(!is.list(glm_options) || (length(glm_options) > 0 && is.null(names(glm_options))) ||
     !all(names(glm_options) %in% c("pMX", "pMXZ"))){
    stop("glm_options must be a list with entries pMX and/or pMXZ, e.g. list(pMX = \".\", pMXZ = \".\").")
  }
  for(nm in names(glm_options)){
    if(!is.null(glm_options[[nm]])){
      if(!is.character(glm_options[[nm]]) || length(glm_options[[nm]]) != 1){
        stop("glm_options$", nm, " must be one character string, the right-hand side of the formula (e.g. \".\").")
      }
      out[[nm]] <- glm_options[[nm]]
    }
  }
  out
}

# check Delta_M_threshold, and warn when Delta_M looks like a cut on M but no threshold is given
check_Delta_M_threshold <- function(Delta_M_threshold, M, Delta_M, uses_gamlss){
  if(is.null(M)){
    return(NULL)
  }
  if(!is.null(Delta_M_threshold)){
    if(!is.numeric(Delta_M_threshold) || length(Delta_M_threshold) != 1 ||
       is.na(Delta_M_threshold) || Delta_M_threshold <= 0){
      stop("Delta_M_threshold must be a single positive number, e.g. 0.2 when Delta_M = as.numeric(M < 0.2).")
    }
    if(!uses_gamlss){
      message("Delta_M_threshold is only used by density_method = \"gamlss\"; it is ignored for GLM and HAL densities.")
      return(NULL)
    }
    if(any(M[Delta_M == 1] >= Delta_M_threshold, na.rm = TRUE)){
      stop("Some participants with Delta_M = 1 have M >= Delta_M_threshold. ",
           "Delta_M_threshold = c assumes Delta_M = 1 only when M < c.")
    }
    if(any(Delta_M == 0 & !is.na(M) & M < Delta_M_threshold)){
      warning("Some participants with Delta_M = 0 have M < Delta_M_threshold, so Delta_M is not only a cut on M. ",
              "The truncated model for pMXD and pMXZD is still fitted on (0, Delta_M_threshold).")
    }
    return(Delta_M_threshold)
  }
  if(uses_gamlss){
    m1 <- M[Delta_M == 1 & !is.na(M)]
    m0 <- M[Delta_M == 0 & !is.na(M)]
    if(length(m1) > 0 && length(m0) > 0 && max(m1) < min(m0)){
      warning(sprintf(paste0(
        "Delta_M looks like a threshold on M: every participant with Delta_M = 1 has M <= %.4g and every participant ",
        "with Delta_M = 0 and observed M has M >= %.4g. If Delta_M = as.numeric(M < c), set Delta_M_threshold = c ",
        "so that pMXD and pMXZD are fitted with a right-truncated GG. An untruncated GG on a sample cut at c ",
        "often fails to converge or gives biased densities."), max(m1), min(m0)), call. = FALSE)
    }
  }
  NULL
}

# ---------------------------------------------------------------------------
# right-truncated GG family
# ---------------------------------------------------------------------------

#' Right-truncated generalized gamma family for GAMLSS
#'
#' @details Builds a \code{gamlss.family} object for the generalized gamma (GG) distribution
#' right-truncated at \code{threshold}, with density
#' \eqn{f(m)\,1(m < c) / F(c)}, where \eqn{f} and \eqn{F} are the GG density and distribution
#' function from \code{gamlss.dist}. MoCo uses it for p(m | a, x, Delta_M = 1) and
#' p(m | a, x, z, Delta_M = 1) when \code{Delta_M = 1(M < threshold)}. The construction follows
#' \code{gamlss.tr::gen.trun(type = "right")}: the score is the GG score minus the derivative of
#' \eqn{\log F(c)}, computed numerically, and the GG expected information is used as the working weight.
#' The threshold is written into the family functions, so the object is self-contained.
#'
#' @param threshold A positive number, the truncation point c.
#' @param mu.link,sigma.link,nu.link Link functions, as in \code{gamlss.dist::GG}.
#'
#' @return A \code{gamlss.family} object that can be passed to \code{gamlss::gamlss(family = )}.
#' @examples
#' if (requireNamespace("gamlss", quietly = TRUE) &&
#'     requireNamespace("gamlss.dist", quietly = TRUE)) {
#'   # motion observed only below a cut-off of 0.2 mm
#'   set.seed(2)
#'   m <- gamlss.dist::rGG(500, mu = 0.12, sigma = 0.5, nu = 0.5)
#'   m <- m[m < 0.2]
#'   fit <- gamlss::gamlss(m ~ 1, family = GGtr(threshold = 0.2),
#'                         control = gamlss::gamlss.control(n.cyc = 200, trace = FALSE))
#'   c(mu = exp(coef(fit, "mu")), sigma = exp(coef(fit, "sigma")), nu = coef(fit, "nu"))
#' }
#' @export
GGtr <- function(threshold, mu.link = "log", sigma.link = "log", nu.link = "identity"){
  if(missing(threshold) || !is.numeric(threshold) || length(threshold) != 1 || threshold <= 0){
    stop("GGtr(): threshold must be a single positive number.")
  }
  if(!requireNamespace("gamlss.dist", quietly = TRUE)){
    stop("Package 'gamlss.dist' is required for GGtr().")
  }
  cc <- threshold
  fam <- do.call(gamlss.dist::GG, list(mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link))

  # log F(c), the truncation constant
  logF <- bquote(function(mu, sigma, nu){
    log(pmax(gamlss.dist::pGG(.(cc), mu = mu, sigma = sigma, nu = nu), 1e-300))
  })
  # central-difference derivative of log F(c) with respect to one parameter
  nd_logF <- function(par){
    args <- list(mu = quote(mu), sigma = quote(sigma), nu = quote(nu))
    plus <- args; plus[[par]] <- bquote(.(as.name(par)) + h)
    minus <- args; minus[[par]] <- bquote(.(as.name(par)) - h)
    bquote({
      h <- 1e-4 * pmax(abs(.(as.name(par))), 1e-2)
      lF <- .(logF)
      (lF(.(plus$mu), .(plus$sigma), .(plus$nu)) - lF(.(minus$mu), .(minus$sigma), .(minus$nu))) / (2 * h)
    })
  }
  make_fun <- function(args, body){
    f <- function() NULL
    formals(f) <- args
    body(f) <- body
    environment(f) <- baseenv()
    f
  }
  a4 <- alist(y = , mu = , sigma = , nu = )

  fam$family <- c("GGtr", paste0("right-truncated Generalized Gamma at ", signif(cc, 6)))
  fam$dldm <- make_fun(a4, bquote(gamlss.dist::GG()$dldm(y, mu, sigma, nu) - .(nd_logF("mu"))))
  fam$dldd <- make_fun(a4, bquote(gamlss.dist::GG()$dldd(y, mu, sigma, nu) - .(nd_logF("sigma"))))
  fam$dldv <- make_fun(a4, bquote(gamlss.dist::GG()$dldv(y, mu, sigma, nu) - .(nd_logF("nu"))))
  # second derivatives: the untruncated GG expected information (kept from gamlss.dist::GG()).
  # The score is exact for the truncated likelihood, so the fixed point is the truncated MLE.
  fam$G.dev.incr <- make_fun(alist(y = , mu = , sigma = , nu = , ... = ),
                             bquote(-2 * (gamlss.dist::dGG(y, mu, sigma, nu, log = TRUE) - .(logF)(mu, sigma, nu))))
  # randomized quantile residuals; rqres() parses pfun as text, so the truncated cdf is passed as text
  pfun_txt <- sprintf(paste0("function(q, mu, sigma, nu) pmin(gamlss.dist::pGG(pmin(q, %.17g), mu = mu, sigma = sigma, nu = nu) / ",
                             "gamlss.dist::pGG(%.17g, mu = mu, sigma = sigma, nu = nu), 1)"), cc, cc)
  fam$rqres <- bquote(rqres(pfun = .(pfun_txt), type = "Continuous", y = y, mu = mu, sigma = sigma, nu = nu))
  fam$y.valid <- make_fun(alist(y = ), bquote(all(y > 0 & y < .(cc))))
  fam$threshold <- cc
  fam
}

# density of the right-truncated GG; 0 at m >= threshold
dGGtr_moco <- function(x, mu, sigma, nu, threshold){
  d <- gamlss.dist::dGG(pmin(x, threshold * (1 - 1e-12)), mu = mu, sigma = sigma, nu = nu) /
    pmax(gamlss.dist::pGG(threshold, mu = mu, sigma = sigma, nu = nu), 1e-300)
  d[x >= threshold] <- 0
  d
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
  user <- if(!is.null(opts$formula)) opts$formula[[label]] else NULL

  as_rhs <- function(f){
    if(inherits(f, "formula")) paste(deparse(f[[length(f)]], width.cutoff = 500L), collapse = " ") else f
  }
  mu_rhs <- if(!is.null(user$mu)) as_rhs(user$mu) else if(inherits(opts$mu, "formula")) as_rhs(opts$mu) else gg_rhs(opts$mu, covars, cont)
  sigma_rhs <- if(!is.null(user$sigma)) as_rhs(user$sigma) else if(inherits(opts$sigma, "formula")) as_rhs(opts$sigma) else gg_rhs(opts$sigma, covars, cont)
  nu_spec <- if(!is.null(user$nu)) user$nu else opts$nu
  nu_fixed <- identical(suppressWarnings(as.numeric(nu_spec)), 0)
  nu_rhs <- if(nu_fixed) "1" else as_rhs(if(identical(nu_spec, "1")) "1" else nu_spec)

  # every variable used must be in W
  used <- unique(c(all.vars(stats::as.formula(paste("~", mu_rhs))),
                   all.vars(stats::as.formula(paste("~", sigma_rhs))),
                   all.vars(stats::as.formula(paste("~", nu_rhs)))))
  missing_vars <- setdiff(used, names(W))
  if(length(missing_vars) > 0){
    stop("gamlss formula for ", label, " uses variables that are not available for this density: ",
         paste(missing_vars, collapse = ", "), ". Note that pMX and pMXD only use A and X.")
  }
  list(mu = paste(".motion ~", mu_rhs), sigma = paste("~", sigma_rhs), nu = paste("~", nu_rhs),
       nu_fixed = nu_fixed, linear_rhs = gg_rhs("linear", covars, cont))
}

# ---------------------------------------------------------------------------
# fitting
# ---------------------------------------------------------------------------

# one call to gamlss(); returns the fit or the error message
gg_try_fit <- function(dat, forms, family, algo, nu_fixed, opts){
  env <- new.env(parent = asNamespace("gamlss"))
  env$.dat <- dat
  mk <- function(txt){ f <- stats::as.formula(txt); environment(f) <- env; f }
  method <- switch(algo,
                   RS = call("RS"),
                   CG = call("CG"),
                   mixed = call("mixed", 10, opts$n.cyc))
  args <- list(formula = mk(forms$mu), sigma.formula = mk(forms$sigma), nu.formula = mk(forms$nu),
               family = family, data = quote(.dat), method = method,
               control = gamlss::gamlss.control(n.cyc = opts$n.cyc, c.crit = opts$c.crit, trace = opts$trace))
  if(nu_fixed){
    args$nu.fix <- TRUE
    args$nu.start <- 0
  }
  cl <- as.call(c(list(as.name("gamlss")), args))
  out <- tryCatch(
    suppressWarnings(eval(cl, env)),
    error = function(e) conditionMessage(e)
  )
  out
}

# fit one density component with the fallback chain:
# each algorithm in opts$method; then nu fixed at 0 (if nu is estimated and nu0_fallback);
# then the GLM log-normal density, with a warning.
gg_component <- function(m, W, label, opts, threshold = NULL){
  W <- gg_clean_data(W)
  dat <- data.frame(.motion = m, W)
  forms <- gg_formulas(W, label, opts)
  family <- if(is.null(threshold)) gamlss.dist::GG() else GGtr(threshold)

  attempts <- lapply(opts$method, function(a) list(algo = a, nu_fixed = forms$nu_fixed))
  if(!forms$nu_fixed && isTRUE(opts$nu0_fallback)){
    attempts <- c(attempts, list(list(algo = "RS", nu_fixed = TRUE), list(algo = "CG", nu_fixed = TRUE)))
  }

  messages <- character(0)
  for(att in attempts){
    fit <- gg_try_fit(dat, forms, family, att$algo, att$nu_fixed, opts)
    if(inherits(fit, "gamlss") && isTRUE(fit$converged) && is.finite(fit$G.deviance)){
      log_row <- data.frame(density = label, method = att$algo, nu_fixed = att$nu_fixed,
                            truncated = !is.null(threshold), converged = TRUE,
                            iterations = fit$iter, nu = mean(fit$nu.fv), n = nrow(dat),
                            fallback = (att$algo != opts$method[1] || att$nu_fixed != forms$nu_fixed),
                            note = paste(messages, collapse = "; "), stringsAsFactors = FALSE)
      return(list(type = "gg", fit = fit, data = dat, threshold = threshold,
                  range_M = range(m), nu_fixed = att$nu_fixed, log = log_row))
    }
    messages <- c(messages, paste0(att$algo, if(att$nu_fixed) "(nu=0)" else "", ": ",
                                   if(is.character(fit)) substr(fit, 1, 80) else "not converged"))
  }

  # last resort: GLM log-normal density on the same rows
  warning("GAMLSS did not converge for ", label, " (", paste(messages, collapse = "; "),
          "). The GLM log-normal density is used for this component.", call. = FALSE)
  lm_fit <- stats::lm(stats::as.formula(paste("log(.motion) ~", forms$linear_rhs)), data = dat)
  log_row <- data.frame(density = label, method = "glm_fallback", nu_fixed = NA,
                        truncated = !is.null(threshold), converged = FALSE, iterations = NA,
                        nu = 0, n = nrow(dat), fallback = TRUE,
                        note = paste(messages, collapse = "; "), stringsAsFactors = FALSE)
  list(type = "glm", fit = lm_fit, data = dat, threshold = threshold,
       range_M = range(m), log = log_row)
}

# evaluate a fitted component at motion values m and covariates W.
# zero_outside = TRUE sets the density to 0 outside the training range of M
# (used for the untruncated pMXD/pMXZD, the same rule as the GLM and HAL densities).
gg_density <- function(obj, m, W, opts, zero_outside = FALSE){
  if(length(m) == 0){
    return(numeric(0))
  }
  W <- gg_clean_data(W)
  newdat <- data.frame(.motion = m, W)
  cc <- obj$threshold

  if(obj$type == "gg"){
    # predict() per parameter; predictAll() looks the family up by name, which a truncated family cannot support
    pars <- lapply(c(mu = "mu", sigma = "sigma", nu = "nu"), function(par){
      if(par == "nu" && isTRUE(obj$nu_fixed)){
        return(rep(0, length(m)))   # nu fixed at 0 has no nu formula to predict from
      }
      suppressWarnings(as.numeric(stats::predict(obj$fit, what = par, newdata = newdat, data = obj$data, type = "response")))
    })
    if(is.null(cc)){
      d <- gamlss.dist::dGG(m, mu = pars$mu, sigma = pars$sigma, nu = pars$nu)
    }else{
      d <- dGGtr_moco(m, mu = pars$mu, sigma = pars$sigma, nu = pars$nu, threshold = cc)
    }
  }else{
    mean_log <- stats::predict(obj$fit, newdata = newdat)
    s <- stats::sd(obj$fit$residuals)
    d <- stats::dnorm(log(m), mean = mean_log, sd = s) / m
    if(!is.null(cc)){
      d <- d / stats::pnorm(log(cc), mean = mean_log, sd = s)
      d[m >= cc] <- 0
    }
  }
  if(zero_outside && is.null(cc)){
    d[m < obj$range_M[1] | m > obj$range_M[2]] <- 0
  }
  d
}

# summary of the density-ratio weights pMXD / pMXZ that enter the efficient influence function
weight_summary <- function(w){
  w <- w[is.finite(w)]
  if(length(w) == 0){
    return(c(n = 0, max_weight = NA, ess_ratio = NA))
  }
  c(n = length(w), max_weight = max(w), ess_ratio = sum(w)^2 / (length(w) * sum(w^2)))
}
