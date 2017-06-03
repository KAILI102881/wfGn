\name{fspecPFGN}
\alias{fspecPFGN}

\title{ Spectral Density of a perturbed Fractional Gaussian Noise }
\description{
Generation of the spectral density of a fractional Gaussian noise contaminated by white noise.
}
\usage{
fspecPFGN(eta, m, SNR=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{ a positive value of the Hurst exponent which is less than 1. }
  \item{m}{ the number of time points. }
  \item{SNR}{ the signal-to-noise ratio. }
}
\details{
This function produces the spectral density of a fractional Gaussian noise contaminated by white noise.
}
\value{
  \item{fspec}{a vector of spectral density with length of the largest integer less than \code{(m-1)/2}.}
  \item{theta1}{a value of the first component of theta.}
}
\references{ 
Wonsang You (2010) Modified Whittle's Maximum Likelihood Estimator for Fractional Gaussian Noises Contaminated by Additive Noises, Technical Reports of the Leibniz Institute for Neurobiology, TR10015.

Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall.
}
\author{ Wonsang You }

\seealso{ \code{\link{fspecFGN}}, \code{\link{fspecPFGN}}, \code{\link{wFGN}} }
\examples{
eta<-0.7; m<-1000; SNR=0
fspec<-fspecPFGN(eta,m,SNR)
}