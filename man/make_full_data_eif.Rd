% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{make_full_data_eif}
\alias{make_full_data_eif}
\title{Evaluate EIF for stochastic interventional (in)direct effects}
\usage{
make_full_data_eif(
  a,
  A,
  Delta_Y,
  Delta_M,
  Y,
  gA,
  gDM,
  gDY_AXZ,
  eta_AXZ,
  eta_AXM,
  xi_AX,
  mu,
  pMXD,
  pMXZ
)
}
\arguments{
\item{a}{Either 0 or 1.}

\item{A}{A binary vector of length n \code{number of participants} representing the diagnosis status.}

\item{Y}{A vector of continuous outcome of interest.}

\item{gA}{Estimate of P(A = 1 | X_i), i = 1,...,n.}

\item{eta_AXZ}{Estimate of E(mu_pseudo_A | A_i, Delta_i, X_i).}

\item{eta_AXM}{Estimate of E(mu_pseudo_A_star | A_i, M_i, X_i).}

\item{xi_AX}{Estimate of E(eta_AXZ | A_i, Delta_i, X_i).}

\item{mu}{Estimate of E(Y| A_i, M_i, X_i, Z_i).}

\item{pMXZ}{Estimate of p(M | A_i, X_i, Z_i).}

\item{Delta}{A binary vector indicating data usability, such as whether data pass a motion criteria.}

\item{gD}{Estimate of P(Delta = 1 | A_i, X_i), i = 1,...,n.}

\item{pMX}{Estimate of p(M | 0, Delta_i = 1, X_i).}
}
\value{
An n-length vector of the estimated EIF evaluated on the observations.
}
\description{
Evaluate EIF for stochastic interventional (in)direct effects
}
