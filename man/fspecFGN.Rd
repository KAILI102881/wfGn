\name{fspecFGN}
\alias{fspecFGN}

\title{ Spectral Density of a Fractional Gaussian Noise }
\description{
Generation of the spectral density of a fractional Gaussian noise.
}
\usage{
fspecFGN(eta,m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{ a positive value of the Hurst exponent which is less than 1. }
  \item{m}{ the number of time points. }
}
\details{
This function produces the spectral density of a fractional Gaussian noise.
}
\value{
  \item{fspec}{a vector of spectral density with length of the largest integer less than \code{(m-1)/2}.}
  \item{theta1}{a value of the first component of theta.}
}
\references{ 
Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall.
}
\author{
Jan Beran (original) and Wonsang You (fine tuning)
}

\seealso{ \code{\link{fspecBFGN}}, \code{\link{fspecPFGN}}, \code{\link{wFGN}} }
\examples{
eta<-0.7; m<-1000
fspec<-fspecFGN(eta,m)
}