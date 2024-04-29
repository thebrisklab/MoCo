#' Simultaneous Confidence Interval Hypothesis Testing
#' 
#' @details Conducts hypothesis testing for significant associations using simultaneous confidence intervals.
#' 
#' @param result The result of fitting based on the one-step function or one-step cross function.
#' @param fwer Family-wise error rates (FWER) for control.
#' @param seed A numeric value to control the randomness of the random variable sampling.
#' 
#' @return A list with named entries:
#'   \item{conf_band}{Simultaneous confidence band based on specified FWER.}
#'   \item{z_score}{A p-length vector of z-scores for each functional connectivity.}
#'   \item{significant_regions}{A p-length vector of TRUE or FALSE, indicating significance for each functional connectivity.}

hypo_test <- function(
    result,
    fwer = 0.05,
    seed = 1
){
  
  # extract information from result
  adj_association <- result$adj_association
  eif_mat <- result$eif_mat
  cov_mat <- result$cov_mat
  
  # z-score calculation
  cov_diff <- sapply(cov_mat, function(mat){
    mat[1,1] - mat[1,2] - mat[2,1] + mat[2,2]
  })
  z_score <- adj_association / sqrt(cov_diff)
  
  # calculate efficient influence function difference
  eif_diff <- do.call(rbind, lapply(eif_mat, function(x){
    x[,2] - x[,1]
  }))
  # number of regions
  p <- nrow(eif_diff)
  # calculate confidence band
  cov_diff_mat <- cor(t(eif_diff))
  # sample from N(0, cov_mat) multivariate normal distribution
  set.seed(seed)
  samples <- MASS::mvrnorm(n = 100000, mu = rep(0, p), Sigma = cov_diff_mat, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  # take max abs value of each sample
  samples_max <- apply(samples, MARGIN = 1, function(x){max(abs(x))})
  
  # number of fwer under consideration
  n_fwer <- length(fwer)
  conf_band <- numeric(n_fwer)
  significant_regions <- matrix(nrow = n_fwer, ncol = length(z_score))
  
  for(i in 1:n_fwer){
    # calculate confidence band
    conf_band[i] <- quantile(samples_max, 1 - fwer[i]) 
    # regions with significant associations
    significant_regions[i,] <- (abs(z_score) > conf_band[i])
  }
  
  return(list(
    z_score = z_score,
    conf_band = conf_band,
    significant_regions = significant_regions)
  )
}




