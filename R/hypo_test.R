#' Simultaneous Confidence Interval Hypothesis Testing
#'
#' @details Both methods simulate the maximum absolute component of a Gaussian
#'   vector whose covariance is the empirical correlation of the EIF contrasts.
#'   \code{"mvn"} draws from a multivariate normal distribution with that
#'   correlation matrix (\code{MASS::mvrnorm}). \code{"multiplier"} uses
#'   independent standard Gaussian observation weights shared across FCs: each
#'   draw is \eqn{\sum_i e_{ij} g_i / \sqrt{n \hat\sigma^2_j}}, where
#'   \eqn{e_{ij}} are the centered EIF contrasts and
#'   \eqn{\hat\sigma^2_j = \sum_i e_{ij}^2 / n}. No step-down or FDX
#'   augmentation is performed. The caller's RNG state is restored. Estimates
#'   and z-statistics are unaffected by the simulation choice; finite Monte
#'   Carlo cutoffs may differ.
#'
#' @param result One-step result with adj_association, eif_mat and cov_mat.
#'   EIF rows must identify the same participants in the same order across FCs.
#' @param fwer FWER probabilities strictly between zero and one.
#' @param method Gaussian simulation method: \code{"multiplier"} (default) or
#'   \code{"mvn"}.
#' @param n_sim Number of simulated Gaussian vectors. Default 10000.
#' @param seed Seed for critical-value simulation, separate from nuisance fitting.
#' @return A list with \code{z_score}, \code{fwer_critical_value}, and
#'   \code{significant_regions}.
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
hypo_test <- function(result,
                      fwer = 0.05,
                      method = c("multiplier", "mvn"),
                      n_sim = 10000L,
                      seed = 1) {

  method <- match.arg(method)

  if (!is.numeric(n_sim) || length(n_sim) != 1L || !is.finite(n_sim) || n_sim < 1 ||
      n_sim != floor(n_sim) || n_sim > .Machine$integer.max)
    stop("n_sim must be one positive integer.", call. = FALSE)
  n_sim <- as.integer(n_sim)

  if (!is.numeric(fwer) || !length(fwer) || any(!is.finite(fwer)) ||
      any(fwer <= 0 | fwer >= 1))
    stop("fwer must contain probabilities strictly between zero and one.", call. = FALSE)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed != floor(seed) || seed > .Machine$integer.max)
    stop("seed must be one nonnegative integer.", call. = FALSE)

  # extract information from result
  adj_association <- result$adj_association
  if (!is.numeric(adj_association) || !length(adj_association) ||
      any(!is.finite(adj_association)))
    stop("adj_association must be a finite numeric vector.", call. = FALSE)

  # number of regions
  p <- length(adj_association)

  if (!is.list(result$eif_mat) || length(result$eif_mat) != p ||
      !is.list(result$cov_mat) || length(result$cov_mat) != p)
    stop("eif_mat and cov_mat must each contain one entry per FC.", call. = FALSE)

  valid_eif <- vapply(result$eif_mat, function(x)
    is.matrix(x) && is.numeric(x) && ncol(x) == 2L && nrow(x) >= 2L && all(is.finite(x)), logical(1))
  if (!all(valid_eif) || length(unique(vapply(result$eif_mat, nrow, integer(1)))) != 1L)
    stop("Every EIF must be a finite n-by-2 matrix with a common n >= 2.", call. = FALSE)

  valid_cov <- vapply(result$cov_mat, function(x)
    is.matrix(x) && is.numeric(x) && identical(dim(x), c(2L, 2L)) && all(is.finite(x)), logical(1))
  if (!all(valid_cov))
    stop("Each covariance must be a finite 2-by-2 matrix.", call. = FALSE)

  # z-score calculation
  cov_diff <- vapply(result$cov_mat, function(x)
    x[1, 1] - x[1, 2] - x[2, 1] + x[2, 2], numeric(1))
  if (any(cov_diff <= 0))
    stop("Every FC contrast must have positive variance.", call. = FALSE)
  z_score <- adj_association / sqrt(cov_diff)

  # efficient influence function difference: p x n (rows = FCs, columns = participants)
  eif_diff <- do.call(rbind, lapply(result$eif_mat, function(x) x[, 2] - x[, 1]))
  n <- ncol(eif_diff)
  # centered EIF contrasts and their variances sum(e^2) / n
  dat_eif <- eif_diff - rowMeans(eif_diff)
  theta_var <- rowSums(dat_eif^2) / n
  if (any(!is.finite(theta_var)) || any(theta_var <= 0))
    stop("Every EIF contrast must have positive finite empirical variance.", call. = FALSE)

  # Local simulation seed: do not overwrite the caller's RNG stream.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  if (method == "mvn") {
    # sample from N(0, cor) multivariate normal distribution
    cov_diff_mat <- stats::cor(t(eif_diff))
    samples <- MASS::mvrnorm(n = n_sim, mu = rep(0, p), Sigma = cov_diff_mat, tol = 1e-6)
    samples <- matrix(samples, nrow = n_sim, ncol = p)
  } else {
    # Weights are generated 1000 draws at a time.
    scaled <- t(dat_eif) / rep(sqrt(n * theta_var), each = n)   # n x p
    samples <- matrix(0, nrow = n_sim, ncol = p)
    first <- 1L
    while (first <= n_sim) {
      count <- min(1000L, n_sim - first + 1L)
      g <- matrix(stats::rnorm(count * n), nrow = count, ncol = n, byrow = TRUE)
      samples[first:(first + count - 1L), ] <- g %*% scaled
      first <- first + count
    }
  }

  # take max abs value of each sample
  samples_max <- apply(abs(samples), 1L, max)

  # calculate confidence band and regions with significant associations
  critical <- as.numeric(stats::quantile(samples_max, probs = 1 - fwer, names = FALSE))
  significant <- vapply(critical, function(cut) abs(z_score) > cut, logical(p))
  significant <- t(matrix(significant, nrow = p, ncol = length(critical)))

  list(z_score = z_score,
       fwer_critical_value = critical,
       significant_regions = significant)
}
