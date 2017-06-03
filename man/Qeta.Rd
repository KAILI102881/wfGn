\name{Qeta}
\alias{Qeta}

\title{ Qeta function }
\description{
Function to be minimized to obtain the maximum likelihood of a fractional Gaussian noise.
}
\usage{
Qeta(eta,n,yper,snr,pertype)
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
Let \eqn{I(f_i)} and \eqn{S(f_i)} be respectively the periodogram of a given perturbed fractional Gaussian noise and the spectral density of perturbed fGn with Hurst exponent \code{eta} and the signal-to-noise ratio \code{SNR} where \eqn{f_i=2\pi *i/n} with \eqn{i=1,...,(n-1)/2}. Then, the value is determined as

\eqn{B=\frac{2\pi}{n}\sum_{i}2\times \left [ \frac{I(f_i)}{S(f_i)} \right ]^2}.

Also, \eqn{A} and \eqn{Tn} are defined as follows.

\eqn{A=\frac{2\pi}{n}\sum_{i}2\times \left [ \frac{I(f_i)}{S(f_i)} \right ]},
\eqn{T_n=\frac{A}{B^2}}

The most parts of this function were adopted from the S-PLUS codes originally developed by Jan Beran. See Beran (1994) for details.
}
\value{
  \item{A,B,Tn}{defined in the above section}
  \item{z}{the test statistics}
  \item{pval}{the p-value}
  \item{fspec}{a vector of spectral density with length of the largest integer less than \code{(m-1)/2}.}
  \item{theta1}{a value of the first component of theta.}
  \item{value}{ a value for minimization }
}
\references{ 
Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall.
}
\author{ Jan Beran (original) and Wonsang You (modifying) }

\seealso{ \code{\link{Csum}}, \code{\link{Lpvar}}, \code{\link{QLCfun}} }
\examples{
n<-1000; H<-0.7; SNR<-10
ts <- perturbFGN(n,H,type="WN",SNR=SNR)
yper<-per(ts)
ts.qeta<-Qeta(H,n,yper,snr=SNR,pertype="per")
}