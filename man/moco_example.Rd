\docType{data}
\name{moco_example}
\alias{moco_example}
\title{Simulated ASD seed-based connectivity data for MoCo}
\format{
A list with 7 elements:
\describe{
\item{A}{Diagnosis, 1 = ASD, 0 = non-ASD.}
\item{X}{Data frame of baseline covariates: \code{X1} sex (1 = male), \code{X2} age (years), \code{X3} handedness (1 = right).}
\item{Z}{Data frame of post-diagnosis covariates: \code{Z1} ADOS total (0 for non-ASD), \code{Z2} full-scale IQ,
\code{Z3} stimulant medication, \code{Z4} non-stimulant medication.}
\item{M}{Mean framewise displacement (mm).}
\item{Delta_M}{Motion usability indicator (1 = usable); not a fixed cut-off on \code{M}.}
\item{Delta_Y}{Imaging-outcome availability (1 = available).}
\item{Y}{400 x 7 matrix of Fisher-z connectivity; rows with \code{Delta_Y = 0} are \code{NA}, and the
seed-to-self column \code{Default} is all \code{NA}.}
}
}
\source{
Simulated. Design motivated by Ran et al. (MoCo) and calibrated to published summaries of
head motion and quality control in children with and without ASD (Nebel et al. 2022, NeuroImage).
}
\usage{
moco_example
}
\description{
A simulated data set with 400 children that mimics a seed-based resting-state
analysis of autism spectrum disorder (ASD), motivated by the Autism Brain
Imaging Data Exchange (ABIDE). The seed is the default mode network, and the
outcomes are its Fisher-z connectivity with six Yeo networks. Head motion is
generated from a generalized gamma distribution whose location and spread depend
on all baseline and post-diagnosis covariates, and motion usability is probabilistic in mean
framewise displacement.
}
\details{
The true MoCo group differences are 0 for outcomes 1-4 and -0.04845 and -0.06815
for outcomes 5 and 6. The generating script is \code{data-raw/moco_example.R} in the
package source (GitHub repository).
}
\examples{
str(moco_example)
}
\keyword{datasets}
