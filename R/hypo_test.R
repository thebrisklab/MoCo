#' Simultaneous Confidence Interval Hypothesis Testing
#'
#' @details Both methods simulate the maximum absolute component of the same
#'   conditional Gaussian vector with covariance equal to the empirical EIF
#'   correlation. "multiplier" uses centered, unit-Euclidean-norm EIF columns
#'   and independent standard Gaussian observation weights shared across FCs.
#' @param result One-step result with adj_association, eif_mat and cov_mat.
#'   EIF rows must refer to the same participants in the same order across FCs.
#' @param fwer Family-wise error rates strictly between zero and one.
#' @param seed Random seed for critical-value simulation (not nuisance fitting).
#' @param method "mvn" (default) or "multiplier". 
#' @param n_sim Number of Gaussian draws. Default 10000.
#' @param chunk_size Draws processed in each memory-bounded batch. Changing
#'   method or chunk size can change a seeded finite Monte Carlo realization.
#' @return List with z_score, fwer_critical_value and significant_regions.
#' @examples
#' # toy efficient influence functions for three outcomes (two groups each)
#' set.seed(1)
#' n <- 100
#' eif <- lapply(1:3, function(j) { a <- rnorm(n); cbind(a, a + rnorm(n)) })
#' toy <- list(adj_association = c(0.05, 0.30, -0.25),
#'             eif_mat = eif,
#'             cov_mat = lapply(eif, function(e) cov(e) / n))
#' hypo_test(toy, fwer = 0.05, n_sim = 10000)
#' @export
hypo_test <- function(result, fwer = 0.05, seed = 1,
                      n_sim = 10000L, method = c("mvn", "multiplier"),
                      chunk_size = 1000L) {
  method <- match.arg(method)
  positive_integer <- function(x, label) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 1 || x != floor(x) || x > .Machine$integer.max)
      stop(label, " must be one positive integer.", call. = FALSE)
    as.integer(x)
  }
  n_sim <- positive_integer(n_sim, "n_sim")
  chunk_size <- min(positive_integer(chunk_size, "chunk_size"), n_sim)
  if (!is.numeric(fwer) || !length(fwer) || any(!is.finite(fwer)) ||
      any(fwer <= 0 | fwer >= 1))
    stop("fwer must contain probabilities strictly between zero and one.", call. = FALSE)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max)
    stop("seed must be one nonnegative integer.", call. = FALSE)
  adj_association <- result$adj_association
  if (!is.numeric(adj_association) || !length(adj_association) ||
      any(!is.finite(adj_association)))
    stop("adj_association must be a finite numeric vector.", call. = FALSE)
  p <- length(adj_association)
  if (!is.list(result$eif_mat) || length(result$eif_mat) != p ||
      !is.list(result$cov_mat) || length(result$cov_mat) != p)
    stop("eif_mat and cov_mat must each contain one entry per FC.", call. = FALSE)
  valid_eif <- vapply(result$eif_mat, function(x)
    is.matrix(x) && is.numeric(x) && ncol(x) == 2L &&
      nrow(x) >= 2L && all(is.finite(x)), logical(1))
  if (!all(valid_eif) || length(unique(vapply(result$eif_mat, nrow, integer(1)))) != 1L)
    stop("Every EIF must be a finite n-by-2 matrix with a common n >= 2.", call. = FALSE)
  valid_cov <- vapply(result$cov_mat, function(x)
    is.matrix(x) && is.numeric(x) && identical(dim(x), c(2L, 2L)) &&
      all(is.finite(x)), logical(1))
  if (!all(valid_cov)) stop("Each covariance must be a finite 2-by-2 matrix.", call. = FALSE)
  cov_diff <- vapply(result$cov_mat, function(x)
    x[1, 1] - x[1, 2] - x[2, 1] + x[2, 2], numeric(1))
  if (any(cov_diff <= 0)) stop("Every FC contrast must have positive variance.", call. = FALSE)
  z_score <- adj_association / sqrt(cov_diff)
  # n x p matrix: preserve participant pairing and all cross-FC correlations.
  n <- nrow(result$eif_mat[[1]])
  eif_diff <- matrix(unlist(lapply(result$eif_mat, function(x)
    x[, 2] - x[, 1]), use.names = FALSE), nrow = n, ncol = p)
  centered <- sweep(eif_diff, 2L, colMeans(eif_diff), "-")
  norms <- sqrt(colSums(centered^2))
  if (any(!is.finite(norms)) || any(norms <= 0))
    stop("Every EIF contrast must have positive finite empirical variance.", call. = FALSE)
  standardized <- sweep(centered, 2L, norms, "/")
  # crossprod(standardized) == cor(eif_diff). Equivalently, the denominator
  # is sqrt((n - 1) * var(EIF)), not sqrt(n * var(EIF)) for R's var().
  if (method == "mvn") {
    correlation <- crossprod(standardized)
    ee <- eigen(correlation, symmetric = TRUE)
    if (any(ee$values < -1e-6 * max(abs(ee$values))))
      stop("EIF correlation matrix is not positive semidefinite.", call. = FALSE)
    factor <- sweep(ee$vectors, 2L, sqrt(pmax(ee$values, 0)), "*")
  }
  # Local simulation seed: do not overwrite the caller's RNG stream.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  maxima <- numeric(n_sim)
  first <- 1L
  while (first <= n_sim) {
    count <- min(chunk_size, n_sim - first + 1L)
    if (method == "mvn") {
      draw <- matrix(stats::rnorm(count * p), nrow = count) %*% t(factor)
    } else {
      draw <- matrix(stats::rnorm(count * n), nrow = count) %*% standardized
    }
    # Avoid apply()'s repeated R callbacks and storing the full B-by-p matrix.
    row_max <- abs(draw[, 1L])
    if (p > 1L) for (j in 2L:p) row_max <- pmax(row_max, abs(draw[, j]))
    maxima[first:(first + count - 1L)] <- row_max
    first <- first + count
  }
  critical <- as.numeric(stats::quantile(maxima, probs = 1 - fwer,
                                       names = FALSE, type = 7))
  significant <- vapply(critical, function(cut) abs(z_score) > cut, logical(p))
  significant <- t(matrix(significant, nrow = p, ncol = length(critical)))
  list(z_score = z_score, fwer_critical_value = critical,
       significant_regions = significant)
}
