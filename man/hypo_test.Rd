\name{hypo_test}
\alias{hypo_test}
\title{Single-step simultaneous Gaussian inference}
\usage{
hypo_test(result, fwer = 0.05, seed = 1,
  n_sim = 10000L, method = c("mvn", "multiplier"), chunk_size = 1000L)
}
\arguments{
\item{result}{One-step result with adj_association, eif_mat and cov_mat. EIF rows must identify the same participants in the same order across FCs.}
\item{fwer}{FWER probabilities strictly between zero and one.}
\item{seed}{Seed for critical-value simulation, separate from nuisance fitting.}
\item{method}{Gaussian simulation method: \code{"mvn"} (default) or \code{"multiplier"}.}
\item{n_sim}{Number of simulated Gaussian vectors. Default 10000.}
\item{chunk_size}{Vectors processed per batch. Changing chunk size can change finite Monte Carlo results.}
}
\value{A list with \code{z_score}, \code{fwer_critical_value}, and \code{significant_regions}.}
\description{Tests multiple contrasts using a single-step maximum absolute z-statistic.}
\details{Both methods target exactly the same conditional Gaussian law with covariance equal to the empirical correlation of centered EIF contrasts. The multiplier method uses Gaussian observation weights shared across FCs and Euclidean-norm standardization, matching R's correlation normalization. No step-down or FDX augmentation is performed. Only simulated maxima are retained. The caller's RNG state is restored. MVN is often faster when the number of FCs is smaller than the sample size. Estimates and z-statistics are unaffected by the simulation choice; finite Monte Carlo cutoffs may differ.}
\references{Du et al. (2025), Causal Inference for Genomic Data with Multiple Heterogeneous Outcomes, Section 5.1. \doi{10.1080/01621459.2025.2468014}}
\examples{
# toy efficient influence functions for three outcomes (two groups each)
set.seed(1)
n <- 100
eif <- lapply(1:3, function(j) { a <- rnorm(n); cbind(a, a + rnorm(n)) })
toy <- list(adj_association = c(0.05, 0.30, -0.25),
            eif_mat = eif,
            cov_mat = lapply(eif, function(e) cov(e) / n))
hypo_test(toy, fwer = 0.05, n_sim = 10000)
}
