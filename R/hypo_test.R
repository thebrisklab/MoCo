#' Simultaneous confidence-band hypothesis tests
#'
#' Computes outcome-specific standardized statistics and simulation-based
#' simultaneous critical values from the correlation of efficient influence-
#' function (EIF) contrasts. The procedure controls the family-wise error rate
#' over all outcomes supplied in `result`.
#'
#' @param result An aggregated one-step result containing
#'   `adj_association`, `eif_mat`, and `cov_mat`. For repeated cross-fitting
#'   seeds, estimates and participant-level EIFs must be averaged and covariance
#'   recomputed before calling this function. [moco()] performs that aggregation
#'   automatically.
#' @param fwer Numeric vector of target family-wise error rates in `(0, 1)`.
#' @param seed Integer seed for the multivariate-normal simulation.
#' @param n_sim Positive integer number of simulation draws. Larger values
#'   reduce Monte Carlo error at greater computational cost.
#' @param chunk_size Positive integer maximum number of draws generated at one
#'   time. Chunking limits memory use when testing many outcomes or using a
#'   large `n_sim`.
#'
#' @details
#' A one-column outcome is handled as a one-dimensional simultaneous test.
#' With multiple columns, a single joint critical value is calculated from the
#' correlation of their participant-level EIF contrasts. Do not run this
#' function separately for each outcome when family-wise inference across all
#' outcomes is desired.
#'
#' @return A list with `z_score`, `conf_band`, `significant_regions`, and
#'   `correlation_matrix`. `significant_regions` has one row per requested FWER
#'   and one column per tested outcome.
#' @export
hypo_test <- function(
    result,
    fwer = 0.05,
    seed = 1,
    n_sim = 100000L,
    chunk_size = 25000L
) {
  required <- c("adj_association", "eif_mat", "cov_mat")
  missing_components <- setdiff(required, names(result))
  if (length(missing_components) > 0L) {
    stop(
      "result is missing: ", paste(missing_components, collapse = ", "),
      "."
    )
  }
  if (!is.numeric(fwer) || length(fwer) == 0L || anyNA(fwer) ||
      any(!is.finite(fwer)) || any(fwer <= 0 | fwer >= 1)) {
    stop("fwer must contain finite numeric values strictly between 0 and 1.")
  }
  n_sim <- .moco_positive_integer(n_sim, "n_sim")
  chunk_size <- .moco_positive_integer(chunk_size, "chunk_size")

  adj_association <- as.numeric(result$adj_association)
  outcome_names <- names(result$adj_association)
  n_outcomes <- length(adj_association)
  if (n_outcomes == 0L || any(!is.finite(adj_association))) {
    stop("adj_association must contain at least one finite estimate.")
  }
  if (!is.list(result$eif_mat) || length(result$eif_mat) != n_outcomes) {
    stop("eif_mat must be a list with one element per outcome.")
  }
  if (!is.list(result$cov_mat) || length(result$cov_mat) != n_outcomes) {
    stop("cov_mat must be a list with one element per outcome.")
  }

  eif_mat <- lapply(seq_len(n_outcomes), function(j) {
    x <- as.matrix(result$eif_mat[[j]])
    if (ncol(x) != 2L || nrow(x) < 2L || any(!is.finite(x))) {
      stop("Each eif_mat element must be a finite n by 2 matrix.")
    }
    x
  })
  n_observations <- vapply(eif_mat, nrow, integer(1))
  if (length(unique(n_observations)) != 1L) {
    stop("All EIF matrices must contain the same participants in the same order.")
  }

  cov_mat <- lapply(result$cov_mat, function(x) {
    x <- as.matrix(x)
    if (!identical(dim(x), c(2L, 2L)) || any(!is.finite(x))) {
      stop("Each cov_mat element must be a finite 2 by 2 matrix.")
    }
    x
  })
  contrast_variance <- vapply(cov_mat, function(mat) {
    mat[1L, 1L] - mat[1L, 2L] - mat[2L, 1L] + mat[2L, 2L]
  }, numeric(1))
  if (any(!is.finite(contrast_variance)) || any(contrast_variance <= 0)) {
    stop("Each estimated contrast variance must be positive and finite.")
  }
  z_score <- adj_association / sqrt(contrast_variance)
  if (!is.null(outcome_names)) names(z_score) <- outcome_names

  eif_contrast <- do.call(
    rbind,
    lapply(eif_mat, function(x) x[, 2L] - x[, 1L])
  )
  if (n_outcomes == 1L) {
    correlation <- matrix(1, 1L, 1L)
  } else {
    correlation <- stats::cor(t(eif_contrast))
    if (any(!is.finite(correlation))) {
      stop(
        "The EIF correlation matrix contains non-finite values; check for ",
        "outcomes with zero-variance EIF contrasts."
      )
    }
    correlation <- (correlation + t(correlation)) / 2
  }
  if (!is.null(outcome_names)) {
    dimnames(correlation) <- list(outcome_names, outcome_names)
  }

  set.seed(seed)
  simulated_max <- numeric(n_sim)
  first <- 1L
  while (first <= n_sim) {
    last <- min(first + chunk_size - 1L, n_sim)
    n_chunk <- last - first + 1L
    draws <- MASS::mvrnorm(
      n = n_chunk,
      mu = rep(0, n_outcomes),
      Sigma = correlation,
      tol = 1e-6,
      empirical = FALSE,
      EISPACK = FALSE
    )
    draws <- matrix(draws, nrow = n_chunk, ncol = n_outcomes)
    simulated_max[first:last] <- apply(abs(draws), 1L, max)
    first <- last + 1L
  }

  conf_band <- as.numeric(stats::quantile(
    simulated_max,
    probs = 1 - fwer,
    names = FALSE
  ))
  fwer_names <- paste0("FWER_", format(fwer, trim = TRUE))
  names(conf_band) <- fwer_names
  significant_regions <- outer(
    conf_band,
    abs(z_score),
    FUN = function(critical_value, statistic) statistic > critical_value
  )
  dimnames(significant_regions) <- list(fwer_names, outcome_names)

  list(
    z_score = z_score,
    conf_band = conf_band,
    significant_regions = significant_regions,
    correlation_matrix = correlation
  )
}

.moco_positive_integer <- function(x, argument) {
  value <- suppressWarnings(as.numeric(x))
  if (
    length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < 1 || value != floor(value) || value > .Machine$integer.max
  ) {
    stop(argument, " must be one positive integer.")
  }
  as.integer(value)
}
