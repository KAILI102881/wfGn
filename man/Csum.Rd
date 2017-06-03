\name{Csum}
\alias{Csum}

\title{ Cumulative Qeta Difference Function }
\description{
Cumulative Qeta Difference Function to be minimized to estimate the parameters of a fractional Gaussian noise.
}
\usage{
Csum(eta,n,yper,snr,pertype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{ a positive value of the Hurst exponent which is less than 1. }
  \item{n}{ the number of time points. }
  \item{yper}{ a vector of periodogram with length of the largest integer less than \code{(n-1)/2} }
  \item{snr}{ the signal-to-noise ratio. }
  \item{pertype}{ the type of periodogram. Possible modes are \code{"per","taper"}. }
}
\details{
Let \eqn{I(f_i)} and \eqn{S(f_i)} be respectively the periodogram of a given perturbed fractional Gaussian noise and the spectral density of perturbed fGn with Hurst exponent \code{eta} and the signal-to-noise ratio \code{SNR} where \eqn{f_i=2\pi *i/n} with \eqn{i=1,...,(n-1)/2}. Then, the value of cumulative Qeta difference function is determined as

\eqn{\eta =\frac{1}{m^*}\sum_{t=1}^{m^*}\left ( I_t-S_t \right )^2 }

where

\eqn{I_t=\sum_{j=1}^{t}I(f_j)/\sum_{j=1}^{m^*}I(f_j)} and
\eqn{S_t=\sum_{j=1}^{t}S_X(f_j)/\sum_{j=1}^{m^*}S_X(f_j)}.

Some parts of this function were adopted from the S-PLUS codes originally developed by Jan Beran. See Beran (1994) for details.
}
\value{
  \item{A,B,Tn}{defined in \link{Qeta}.}
  \item{z}{the test statistics}
  \item{pval}{the p-value}
  \item{fspec}{a vector of spectral density with length of the largest integer less than \code{(m-1)/2}.}
  \item{theta1}{a value of the first component of theta.}
  \item{value}{ a value for minimization }
}
\references{ 
Wonsang You (2010) Modified Whittle's Maximum Likelihood Estimator for Fractional Gaussian Noises Contaminated by Additive Noises, Technical Reports of the Leibniz Institute for Neurobiology, TR10015.

Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall.
}
\author{ Wonsang You }

\seealso{ \code{\link{Qeta}}, \code{\link{Lpvar}}, \code{\link{QLCfun}} }
\examples{
n<-1000; H<-0.7; SNR<-10
ts <- perturbFGN(n,H,type="WN",SNR=SNR)
yper<-per(ts)
ts.csum<-Csum(H,n,yper,snr=SNR,pertype="per")
}