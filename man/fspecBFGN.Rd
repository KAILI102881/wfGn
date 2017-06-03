\name{fspecBFGN}
\alias{fspecBFGN}

\title{ Spectral Density of a bivariate Fractional Gaussian Noise }
\description{
Generation of the spectral density of a bivariate fractional Gaussian noise.
}
\usage{
fspecBFGN(eta1, eta2, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta1}{ the Hurst exponent of the first component of a bivariate fractional Gaussian noise. }
  \item{eta2}{ the Hurst exponent of the second component of a bivariate fractional Gaussian noise. }
  \item{m}{ the number of time points. }
}
\details{
This function produces the spectral density of a bivariate fractional Gaussian noise.
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
eta1<-0.7; eta2<-0.8; m<-1000
fspec<-fspecBFGN(eta1,eta2,m)
}