\name{CetaFGN}
\alias{CetaFGN}

\title{ Covariance Matrix of fGn Parameters }
\description{
Covariance matrix of parameters \eqn{\widehat{\eta }} of fractional Gaussian noise.
}
\usage{
CetaFGN(eta,snr=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{ a positive value of the Hurst exponent which is less than 1. }
  \item{snr}{ the signal-to-noise ratio. Default is \code{NULL}.}
}
\details{
The most parts of this function were adopted from the S-PLUS codes originally developed by Jan Beran. See Beran (1994) for details.
}
\value{
A covariance matrix of the estimated fGn parameters with size \eqn{M \times M} where \eqn{M} is the length of \eqn{\widehat{\eta }}.
}
\references{ 
Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall.
}
\author{ 
Jan Beran (original) and Wonsang You (fine tuning)
}

\seealso{ \link{fspecFGN} }
\examples{
eta <- c(H=0.7); SNR <- 10
cov <- CetaFGN(eta,SNR)
}