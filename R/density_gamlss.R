# GAMLSS conditional-density helpers for the four motion densities used by MoCo:
#   pMX   = p(M | A, X)                 fit on Delta_Y = 1, evaluated at the observed A on Delta_Y = 1 rows
#   pMXD  = p(M | Delta_M = 1, A, X)    fit on Delta_M = 1, evaluated at A = 0 on every row with observed M (q)
#   pMXZ  = p(M | A, X, Z)              fit on Delta_Y = 1, evaluated at the observed A, A = 0 and A = 1 on Delta_Y = 1 rows
#   pMXZD = p(M | Delta_M = 1, A, X, Z) fit on Delta_M = 1, evaluated at the observed A on every row with observed M
# No density is set to 0 where Delta_M = 0, whether Delta_M was supplied or created from thresh
# (the same evaluation rule as the GLM branch).
# Each density is fitted once per call (once per fold when cross-fitting); all evaluation rows
# (training and validation fold, each A value needed) are stacked in a single predictAll() call.
#
# gamlss_options (all entries optional):
#   family   family name ("GG"), family function (GG) or family object (GG()); or a named list
#            by role with an optional `default`, e.g. list(default = "GG", pMXD = "BCTo").
#            The family must be continuous on (0, Inf). Default "GG".
#   formulas list(mu, sigma, nu, tau). Each entry is a formula, a character string, or a keyword:
#            "linear" = all available columns (A, X and, for pMXZ/pMXZD, Z) entering linearly;
#            "pb"     = the same columns, with pb() on continuous columns (numeric with at least
#                       `continuous_min_unique` distinct values) and binary/factor columns linear;
#            "1"      = intercept only.
#            A "." in a formula expands to all available columns. A per-role list can be given
#            as formulas = list(mu = "pb", sigma = "pb", pMXD = list(sigma = "1")) or with
#            `default`. Defaults: mu = "linear", sigma = nu = tau = "1".
#   method   fitting algorithm passed to gamlss(): "RS()" (default), "CG()", "mixed(1, 50)", as a
#            string or a call such as quote(mixed(1, 50)).
#   control  list of gamlss.control() arguments (n.cyc, c.crit, trace, mu.step, ...) or a
#            gamlss.control object; n.cyc defaults to 200 and trace to FALSE.
#   i.control list of glim.control() arguments (cc, cyc, bf.cyc, ...) or a glim.control object.
#   weights  prior weights: a scalar, or a vector of length n (full data; subset to each fold's fit rows).
#   mu.start, sigma.start, nu.start, tau.start  starting values: scalars, or vectors of length n.
#   mu.fix, sigma.fix, nu.fix, tau.fix  logical, keep that parameter fixed at its start value.
#   args     list of any further arguments forwarded to gamlss::gamlss() (e.g. list(contrasts = ...)).
#   start    "auto" (default): Box-Cox families (BCCGo, BCPEo, BCTo, ...) start from the LOGNO fit with
#            the same mu/sigma formulas (mu = exp(LOGNO mu), sigma = LOGNO sigma, nu = 0.1, tau = 2)
#            unless start values are supplied; "default": gamlss starting values.
#   retry    TRUE (default): if the fit fails or does not converge, refit (1) starting from the fit
#            with the same mu formula and intercept-only sigma/nu/tau, then (2) with CG(), then (3) with
#            RS() if another method was requested. The method
#            used is recorded in gamlss_info. FALSE: no retry.
#   clamp    TRUE (default): continuous covariates of the rows being predicted are limited to the
#            range of the training rows of that density (no extrapolation of pb()/sigma terms);
#            the number of clamped values is reported in gamlss_info$n_clamped. FALSE: no clamping.
#   on_nonconvergence "stop" (default) or "warn" (after any retries).
#   continuous_min_unique  default 10.
#   pMX, pMXD, pMXZ, pMXZD  per-role overrides: a list with any of the entries above (except
#            continuous_min_unique), e.g. pMXD = list(family = "BCTo", method = "CG()", control = list(n.cyc = 400)).
#            control, i.control and args are merged with the global ones; formulas are merged by parameter.

`%||%` <- function(x, y) if (is.null(x)) y else x

.moco_gamlss_roles <- c("pMX", "pMXD", "pMXZ", "pMXZD")
.moco_gamlss_params <- c("mu", "sigma", "nu", "tau")
.moco_gamlss_boxcox <- c("BCCG", "BCCGo", "BCPE", "BCPEo", "BCT", "BCTo")
.moco_gamlss_rowvec <- c("weights", "mu.start", "sigma.start", "nu.start", "tau.start")
.moco_gamlss_fixes <- c("mu.fix", "sigma.fix", "nu.fix", "tau.fix")

# Resolve a family into its object, name, parameters and density function, and check that it
# is a continuous distribution on the positive real line.
.moco_gamlss_family <- function(family) {
  if (inherits(family, "gamlss.family")) {
    obj <- family
  } else {
    fun <- family
    if (is.character(family)) {
      fun <- get0(family, envir = asNamespace("gamlss.dist"), mode = "function") %||%
        get0(family, envir = globalenv(), mode = "function")
      if (is.null(fun)) stop(sprintf("gamlss family '%s' was not found in gamlss.dist or the global environment.", family), call. = FALSE)
    }
    if (!is.function(fun)) stop("gamlss_options$family must be a family name, a family function or a family object.", call. = FALSE)
    obj <- fun()
  }
  name <- obj$family[1]
  dfun <- get0(paste0("d", name), envir = asNamespace("gamlss.dist"), mode = "function") %||%
    get0(paste0("d", name), envir = globalenv(), mode = "function")
  if (is.null(dfun)) stop(sprintf("No density function d%s() found for gamlss family '%s'.", name, name), call. = FALSE)
  yv <- function(y) isTRUE(tryCatch(obj$y.valid(y), error = function(e) FALSE))
  positive <- identical(obj$type, "Continuous") && yv(c(1e-6, 0.5, 1, 1e6)) && !yv(0) && !yv(-1)
  if (!positive) {
    stop(sprintf(paste0("gamlss family '%s' is not a continuous distribution on the positive real line (0, Inf). ",
                        "Motion M needs such a family, e.g. LOGNO, GA, IG, WEI3, GG, BCCGo, BCPEo, BCTo or GB2."), name),
         call. = FALSE)
  }
  list(obj = obj, name = name, params = names(obj$parameters), dfun = dfun)
}

# Columns that enter keyword formulas, and which of them are continuous.
.moco_gamlss_columns <- function(frame, min_unique) {
  vars <- setdiff(names(frame), "M")
  vars <- vars[vapply(frame[vars], function(v) length(unique(v[!is.na(v)])) > 1L, logical(1))]
  cont <- vars[vapply(frame[vars], function(v) is.numeric(v) && length(unique(v[!is.na(v)])) >= min_unique, logical(1))]
  list(vars = vars, cont = cont)
}

.moco_gamlss_name <- function(x) ifelse(make.names(x) == x, x, paste0("`", gsub("`", "\\\\`", x), "`"))

# Make the gamlss additive-term constructors visible to a formula without masking user objects.
.moco_gamlss_formula_env <- function(parent) {
  env <- new.env(parent = parent)
  ns <- asNamespace("gamlss")
  for (nm in c("pb", "pbm", "pbc", "pbz", "pbp", "cs", "scs", "fp", "pp", "random", "re", "ri", "lo", "ps", "pvc")) {
    if (exists(nm, envir = ns, inherits = FALSE) && !exists(nm, envir = parent, inherits = TRUE)) {
      assign(nm, get(nm, envir = ns), envir = env)
    }
  }
  env
}

# Build the formula of one distribution parameter.
.moco_gamlss_formula <- function(f, param, cols) {
  lhs <- if (param == "mu") "M ~ " else "~ "
  parent <- globalenv()
  if (is.null(f)) f <- if (param == "mu") "linear" else "1"
  if (is.character(f) && length(f) == 1L && tolower(trimws(f)) %in% c("linear", "pb", "1", "intercept", "~1")) {
    key <- tolower(trimws(f))
    rhs <- if (key %in% c("1", "intercept", "~1") || !length(cols$vars)) {
      "1"
    } else {
      lin <- setdiff(cols$vars, if (key == "pb") cols$cont)
      sm <- if (key == "pb") sprintf("pb(%s)", .moco_gamlss_name(intersect(cols$vars, cols$cont)))
      paste(c(.moco_gamlss_name(lin), sm), collapse = " + ")
    }
    f <- stats::as.formula(paste0(lhs, rhs), env = parent)
  } else {
    if (inherits(f, "formula")) {
      parent <- environment(f) %||% globalenv()
      txt <- paste(deparse(f, width.cutoff = 500L), collapse = " ")
    } else if (is.character(f) && length(f) == 1L) {
      txt <- f
    } else {
      stop(sprintf("gamlss formula for %s must be a formula, a character string or one of \"linear\", \"pb\", \"1\".", param), call. = FALSE)
    }
    rhs <- if (grepl("~", txt, fixed = TRUE)) sub("^[^~]*~", "", txt) else txt
    if (grepl("(^|[^A-Za-z0-9_.])\\.([^A-Za-z0-9_.]|$)", rhs)) {
      tmpl <- as.data.frame(stats::setNames(rep(list(0), length(cols$vars) + 1L), c("M", cols$vars)), check.names = FALSE)
      expanded <- stats::terms(stats::as.formula(paste("M ~", rhs)), data = tmpl)
      labs <- attr(expanded, "term.labels")
      rhs <- paste(c(if (attr(expanded, "intercept") == 0) "-1", if (length(labs)) labs else "1"), collapse = " + ")
    }
    f <- stats::as.formula(paste0(lhs, rhs), env = parent)
  }
  environment(f) <- .moco_gamlss_formula_env(environment(f))
  f
}

# Options of one role: global options with the role's override list (gamlss_options[[role]]) applied.
.moco_gamlss_role_options <- function(options, role) {
  base <- options[setdiff(names(options), .moco_gamlss_roles)]
  ro <- options[[role]]
  if (is.null(ro)) return(base)
  if (!is.list(ro)) stop(sprintf("gamlss_options$%s must be a list of options for that density.", role), call. = FALSE)
  for (k in names(ro)) {
    if (k %in% c("control", "i.control", "args") && is.list(base[[k]]) && is.list(ro[[k]])) {
      base[[k]] <- utils::modifyList(as.list(unclass(base[[k]])), as.list(unclass(ro[[k]])))
    } else if (k %in% c("formulas", "formula")) {
      fl <- ro[[k]]
      if (inherits(fl, "formula") || is.character(fl)) fl <- list(mu = fl)
      base$.role_formulas <- fl
    } else {
      base[k] <- list(ro[[k]])
    }
  }
  base
}

# Family and formulas of one role.
.moco_gamlss_role_spec <- function(options, role, cols) {
  fam <- options$family
  if (is.list(fam) && !inherits(fam, "gamlss.family")) fam <- fam[[role]] %||% fam$default
  fam <- .moco_gamlss_family(fam %||% "GG")
  fl <- options$formulas %||% options$formula %||% list()
  if (inherits(fl, "formula") || is.character(fl)) fl <- list(mu = fl)
  spec <- fl[intersect(names(fl), .moco_gamlss_params)]
  if (is.list(fl$default)) spec[names(fl$default)] <- fl$default
  if (is.list(fl[[role]])) spec[names(fl[[role]])] <- fl[[role]]
  rf <- options$.role_formulas
  if (is.list(rf)) spec[intersect(names(rf), .moco_gamlss_params)] <- rf[intersect(names(rf), .moco_gamlss_params)]
  forms <- lapply(stats::setNames(.moco_gamlss_params, .moco_gamlss_params),
                  function(p) .moco_gamlss_formula(spec[[p]], p, cols))
  list(family = fam, formulas = forms[intersect(.moco_gamlss_params, fam$params)])
}

# Fill the defaults of the fitting options of one role and check them.
.moco_gamlss_role_fill <- function(ro, role, n) {
  ro$on_nonconvergence <- match.arg(ro$on_nonconvergence %||% "stop", c("stop", "warn"))
  ro$retry <- ro$retry %||% TRUE
  ro$clamp <- ro$clamp %||% TRUE
  ro$start <- match.arg(ro$start %||% "auto", c("auto", "default"))
  ctrl <- utils::modifyList(list(n.cyc = 200, trace = FALSE), as.list(unclass(ro$control %||% list())))
  ro$control <- do.call(gamlss::gamlss.control, ctrl)
  if (!is.null(ro$i.control)) ro$i.control <- do.call(gamlss::glim.control, as.list(unclass(ro$i.control)))
  if (is.character(ro$method)) ro$method <- str2lang(ro$method)
  if (!is.null(ro$args) && (!is.list(ro$args) || is.null(names(ro$args)) || any(names(ro$args) == ""))) {
    stop(sprintf("gamlss_options$args (%s) must be a named list of arguments for gamlss::gamlss().", role), call. = FALSE)
  }
  for (k in .moco_gamlss_rowvec) {
    v <- ro[[k]]
    if (is.null(v)) next
    if (!is.numeric(v) || !(length(v) %in% c(1L, n)) || any(!is.finite(v))) {
      stop(sprintf("gamlss_options$%s (%s) must be a finite numeric scalar or a vector of length n = %d.", k, role, n), call. = FALSE)
    }
    if (k == "weights" && any(v < 0)) stop("gamlss_options$weights must be non-negative.", call. = FALSE)
  }
  for (k in .moco_gamlss_fixes) {
    if (!is.null(ro[[k]]) && !(is.logical(ro[[k]]) && length(ro[[k]]) == 1L && !is.na(ro[[k]]))) {
      stop(sprintf("gamlss_options$%s (%s) must be TRUE or FALSE.", k, role), call. = FALSE)
    }
  }
  ro
}

# Check packages and data once, resolve every role's options, and fix the continuous columns
# from the full data.
.moco_gamlss_prepare <- function(options, X, Z, M) {
  if (isTRUE(options$.prepared)) return(options)
  if (!requireNamespace("gamlss", quietly = TRUE) || !requireNamespace("gamlss.dist", quietly = TRUE)) {
    stop("Install both 'gamlss' and 'gamlss.dist' to use gamlss_options.", call. = FALSE)
  }
  if (!is.list(options)) stop("gamlss_options must be a list.", call. = FALSE)
  M <- as.numeric(unlist(M))
  if (any(!is.na(M) & (!is.finite(M) | M <= 0))) {
    stop("GAMLSS motion densities need every observed motion value to be finite and > 0 (M > 0).", call. = FALSE)
  }
  n <- length(M)
  min_unique <- options$continuous_min_unique %||% 10
  frame_x <- .moco_density_frame(M, rep(0:1, length.out = n), X)
  frame_z <- .moco_density_frame(M, rep(0:1, length.out = n), X, Z)
  cols <- list(x = .moco_gamlss_columns(frame_x, min_unique), z = .moco_gamlss_columns(frame_z, min_unique))
  roles <- list()
  for (role in .moco_gamlss_roles) {
    ro <- .moco_gamlss_role_fill(.moco_gamlss_role_options(options, role), role, n)
    ro$.cols <- if (role %in% c("pMXZ", "pMXZD")) cols$z else cols$x
    ro$.spec <- .moco_gamlss_role_spec(ro, role, ro$.cols)   # validates families and formulas early
    roles[[role]] <- ro
  }
  list(.prepared = TRUE, .n = n, .roles = roles)
}

# Restrict full-data row vectors (weights, *.start of length n) to the rows of one fold.
.moco_gamlss_subset <- function(options, rows) {
  for (role in names(options$.roles)) {
    for (k in .moco_gamlss_rowvec) {
      v <- options$.roles[[role]][[k]]
      if (length(v) > 1L) options$.roles[[role]][[k]] <- v[rows]
    }
  }
  options$.n <- length(rows)
  options
}

# Error message when a GAMLSS density fails or does not converge: names the density, family and
# fold, and suggests other families for that density.
.moco_gamlss_fail_msg <- function(role, family, label, n_fit, what, log) {
  rows <- if (role %in% c("pMXD", "pMXZD")) "Delta_M = 1" else "Delta_Y = 1"
  others <- setdiff(c("LOGNO", "GA", "WEI3", "BCCGo"), family)
  paste0(
    sprintf("GAMLSS density %s (family %s, %s) %s. It was fitted on the %d rows with %s.\n", role, family, label, what, n_fit, rows),
    sprintf("Try another family for this density, e.g. gamlss_options = list(family = list(default = \"%s\", %s = \"%s\")).\n",
            family, role, others[1]),
    sprintf("Other families to try: %s. LOGNO (lognormal) is the most stable and matches the logGLM method; ", paste(others, collapse = ", ")),
    "simpler formulas (e.g. sigma = \"1\") or a larger control$n.cyc can also help.\n",
    "Attempts: ", paste(log, collapse = "; "))
}

.moco_gamlss_fit <- function(train, role, spec, ro, fit_rows, label) {
  args <- list(formula = spec$formulas$mu, sigma.formula = spec$formulas$sigma,
               family = spec$family$obj, data = train, control = ro$control)
  if ("nu" %in% spec$family$params) args$nu.formula <- spec$formulas$nu
  if ("tau" %in% spec$family$params) args$tau.formula <- spec$formulas$tau
  if (!is.null(ro$method)) args$method <- ro$method
  if (!is.null(ro$i.control)) args$i.control <- ro$i.control
  for (k in .moco_gamlss_rowvec) {
    v <- ro[[k]]
    if (!is.null(v)) args[[k]] <- if (length(v) > 1L) v[fit_rows] else v
  }
  if (length(args$weights) == 1L) args$weights <- rep(args$weights, nrow(train))
  for (k in .moco_gamlss_fixes) if (!is.null(ro[[k]])) args[[k]] <- ro[[k]]
  if (length(ro$args)) args[names(ro$args)] <- ro$args
  warn <- character()
  start_note <- NULL
  user_start <- any(paste0(.moco_gamlss_params, ".start") %in% names(args)) || !is.null(args$start.from)
  if (identical(ro$start, "auto") && spec$family$name %in% .moco_gamlss_boxcox && !user_start) {
    # Box-Cox families: start mu and sigma from the LOGNO fit with the same mu/sigma formulas,
    # nu = 0.1 (close to the log transform) and tau = 2 (normal tails)
    la <- args[setdiff(names(args), c("nu.formula", "tau.formula", .moco_gamlss_fixes, "method"))]
    la$family <- gamlss.dist::LOGNO()
    ln <- tryCatch(suppressWarnings(do.call(gamlss::gamlss, la)), error = function(e) NULL)
    if (!is.null(ln)) {
      args$mu.start <- exp(stats::fitted(ln, "mu")); args$sigma.start <- stats::fitted(ln, "sigma")
      if ("nu" %in% spec$family$params) args$nu.start <- 0.1
      if ("tau" %in% spec$family$params) args$tau.start <- 2
      start_note <- "started from the LOGNO fit"
    }
  }
  try_fit <- function(a) withCallingHandlers(
    tryCatch(do.call(gamlss::gamlss, a), error = function(e) e),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  ok <- function(f) !inherits(f, "error") && isTRUE(f$converged)
  t0 <- proc.time()[["elapsed"]]
  fit <- try_fit(args); method <- paste(c(deparse(args$method %||% quote(RS())), start_note), collapse = " ")
  log <- paste0(method, ": ", if (inherits(fit, "error")) conditionMessage(fit) else if (!ok(fit)) "not converged" else "ok")
  if (!ok(fit) && isTRUE(ro$retry)) {
    # retry 1: start from the fit with the same mu formula and intercept-only sigma/nu/tau
    simple <- args
    one <- .moco_gamlss_formula("1", "sigma", NULL)
    simple$sigma.formula <- one
    if (!is.null(simple$nu.formula)) simple$nu.formula <- one
    if (!is.null(simple$tau.formula)) simple$tau.formula <- one
    f1 <- try_fit(simple)
    if (!inherits(f1, "error")) {
      a2 <- args; a2$start.from <- f1
      a2[paste0(.moco_gamlss_params, ".start")] <- NULL
      fit2 <- try_fit(a2)
      log <- c(log, paste0("start from sigma ~ 1 fit: ", if (inherits(fit2, "error")) conditionMessage(fit2) else if (!ok(fit2)) "not converged" else "ok"))
      if (ok(fit2) || inherits(fit, "error")) { fit <- fit2; method <- paste(method, "started from the sigma ~ 1 fit") }
    } else log <- c(log, paste0("sigma ~ 1 start fit: ", conditionMessage(f1)))
    # retry 2: CG algorithm
    if (!ok(fit)) {
      a3 <- args; a3$method <- call("CG", n.cyc = ro$control$n.cyc)
      fit3 <- try_fit(a3)
      log <- c(log, paste0("CG(): ", if (inherits(fit3, "error")) conditionMessage(fit3) else if (!ok(fit3)) "not converged" else "ok"))
      if (ok(fit3) || inherits(fit, "error")) { fit <- fit3; method <- "CG()" }
    }
    # retry 3: RS algorithm, when another method (e.g. mixed() or CG()) was requested
    # (on the example data mixed(1, 50) and CG() fail on the Delta_M = 1 rows, RS() converges)
    if (!ok(fit) && !is.null(args$method) && !grepl("^RS\\(", deparse(args$method)[1])) {
      a4 <- args; a4$method <- NULL
      fit4 <- try_fit(a4)
      log <- c(log, paste0("RS(): ", if (inherits(fit4, "error")) conditionMessage(fit4) else if (!ok(fit4)) "not converged" else "ok"))
      if (ok(fit4) || inherits(fit, "error")) { fit <- fit4; method <- paste(c("RS()", start_note), collapse = " ") }
    }
  }
  if (inherits(fit, "error")) {
    stop(.moco_gamlss_fail_msg(role, spec$family$name, label, nrow(train), "failed to fit", log), call. = FALSE)
  }
  attr(fit, "moco_time") <- proc.time()[["elapsed"]] - t0
  attr(fit, "moco_warnings") <- unique(warn)
  attr(fit, "moco_method") <- method
  attr(fit, "moco_attempts") <- paste(log, collapse = "; ")
  if (!isTRUE(fit$converged)) {
    msg <- .moco_gamlss_fail_msg(role, spec$family$name, label, nrow(train),
                                 sprintf("did not converge in %d cycles", fit$iter), log)
    if (ro$on_nonconvergence == "stop") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  fit
}

# Densities of all rows in `newdata` (one predictAll call).
.moco_gamlss_eval <- function(fit, spec, newdata, train, role, label) {
  if (!nrow(newdata)) return(list(density = numeric(), warn = character()))
  warn <- character()
  pred <- withCallingHandlers(
    gamlss::predictAll(fit, newdata = newdata, data = train, type = "response"),
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  pars <- lapply(spec$family$params, function(p) as.numeric(pred[[p]]))
  names(pars) <- spec$family$params
  density <- do.call(spec$family$dfun, c(list(x = newdata$M), pars))
  if (any(!is.finite(density)) || any(density < 0)) {
    stop(sprintf("Invalid GAMLSS density from %s (%s, %s).", role, spec$family$name, label), call. = FALSE)
  }
  list(density = density, warn = unique(warn))
}

.moco_gamlss_info_row <- function(fit, spec, role, label, pred_warn, time_pred) {
  fw <- attr(fit, "moco_warnings")
  f <- function(p) if (is.null(spec$formulas[[p]])) NA_character_ else paste(deparse(spec$formulas[[p]], width.cutoff = 500L), collapse = " ")
  data.frame(fold = label, role = role, family = spec$family$name,
             mu = f("mu"), sigma = f("sigma"), nu = f("nu"), tau = f("tau"),
             n = fit$N, edf = fit$df.fit, G_deviance = fit$G.deviance,
             AIC = fit$G.deviance + 2 * fit$df.fit, SBC = fit$G.deviance + log(fit$N) * fit$df.fit,
             converged = isTRUE(fit$converged), iter = fit$iter, method = attr(fit, "moco_method"),
             attempts = attr(fit, "moco_attempts"),
             time_fit = attr(fit, "moco_time"), time_pred = time_pred,
             n_warnings = length(fw) + length(pred_warn),
             warnings = paste(substr(unique(c(fw, pred_warn)), 1, 150), collapse = " | "),
             stringsAsFactors = FALSE)
}

# All four motion densities on the training rows and (optionally) the validation rows.
# `train` / `valid` are lists with M, A, X, Z, Delta_Y, Delta_M (valid = NULL: no cross-fitting).
# pMX/pMXZ are fitted and evaluated on Delta_Y = 1 rows; pMXD/pMXZD are fitted on Delta_M = 1 rows and
# evaluated on every row with observed M (pMXD at A = 0, pMXZD at the observed A); nothing is zeroed.
.moco_gamlss_densities <- function(train, valid = NULL, options, label = "full data") {
  sets <- list(train = train, valid = valid)
  sets <- sets[!vapply(sets, is.null, logical(1))]
  frames <- lapply(sets, function(s) list(x = .moco_density_frame(s$M, s$A, s$X), z = .moco_density_frame(s$M, s$A, s$X, s$Z),
                                          DY = as.numeric(unlist(s$Delta_Y)), DM = as.numeric(unlist(s$Delta_M)),
                                          M = as.numeric(unlist(s$M))))
  if (length(frames$train$M) != options$.n) stop("Internal error: GAMLSS options were not subset to the training rows.", call. = FALSE)
  out <- list(); info <- list(); fits <- list()
  roles <- list(pMX = list(z = FALSE, fit_on = "DY", A = "obs"),
                pMXD = list(z = FALSE, fit_on = "DM", A = "0"),
                pMXZ = list(z = TRUE, fit_on = "DY", A = c("obs", "0", "1")),
                pMXZD = list(z = TRUE, fit_on = "DM", A = "obs"))
  for (role in names(roles)) {
    r <- roles[[role]]
    ro <- options$.roles[[role]]
    cols <- ro$.cols
    spec <- ro$.spec
    ftr <- frames$train
    fit_rows <- which(ftr[[r$fit_on]] == 1 & !is.na(ftr$M))
    if (!length(fit_rows)) stop(sprintf("No training rows for GAMLSS %s (%s).", role, label), call. = FALSE)
    fr <- if (r$z) ftr$z else ftr$x
    keep <- c("M", cols$vars)
    trn <- fr[fit_rows, keep, drop = FALSE]
    fit <- .moco_gamlss_fit(trn, role, spec, ro, fit_rows, label)
    # evaluation rows per set
    blocks <- list(); where <- list()
    for (s in names(frames)) {
      fs <- frames[[s]]
      ev <- if (role %in% c("pMX", "pMXZ")) which(fs$DY == 1 & !is.na(fs$M)) else which(!is.na(fs$M))
      base <- (if (r$z) fs$z else fs$x)[ev, keep, drop = FALSE]
      for (a in r$A) {
        b <- base
        if (a != "obs") b$A <- rep(as.numeric(a), nrow(b))
        blocks[[length(blocks) + 1L]] <- b
        where[[length(where) + 1L]] <- list(set = s, a = a, rows = ev)
      }
    }
    newdata <- do.call(rbind, blocks)
    n_clamped <- 0L
    if (isTRUE(ro$clamp)) {
      # keep continuous covariates of the prediction rows inside the training-fold range, so that
      # smooth terms (pb) and sigma models are not extrapolated (ABCD: one T-score of 147 vs a
      # training maximum of 120 gave sigma below its training range and a density ratio of 3e5)
      for (v in intersect(cols$cont, names(newdata))) {
        rg <- range(trn[[v]])
        outside <- newdata[[v]] < rg[1] | newdata[[v]] > rg[2]
        n_clamped <- n_clamped + sum(outside)
        newdata[[v]] <- pmin(pmax(newdata[[v]], rg[1]), rg[2])
      }
    }
    t0 <- proc.time()[["elapsed"]]
    ev <- .moco_gamlss_eval(fit, spec, newdata, trn, role, label)
    tp <- proc.time()[["elapsed"]] - t0
    pos <- 0L
    for (w in where) {
      nm <- paste0(role, "n_", if (w$a == "obs") "A" else paste0("A", w$a), "_", w$set)
      v <- rep(NA_real_, length(frames[[w$set]]$M))
      if (length(w$rows)) v[w$rows] <- ev$density[pos + seq_along(w$rows)]
      pos <- pos + length(w$rows)
      out[[nm]] <- v
    }
    info[[role]] <- .moco_gamlss_info_row(fit, spec, role, label, ev$warn, tp)
    info[[role]]$n_clamped <- n_clamped
    fits[[role]] <- fit
  }
  out$info <- do.call(rbind, info); rownames(out$info) <- NULL
  out$fits <- fits
  out
}

.moco_density_frame <- function(M, A, X, Z = NULL) {
  out <- data.frame(M = as.numeric(unlist(M)), A = as.numeric(unlist(A)), as.data.frame(X), check.names = FALSE)
  if (!is.null(Z)) out <- cbind(out, as.data.frame(Z))
  out
}
