\name{Lpvar}
\alias{Lpvar}

\title{ Log Periodogram Variance Function }
\description{
Log Periodogram Variance Function to be minimized to estimate the parameters of a fractional Gaussian noise.
}
\usage{
Lpvar(eta,n,yper,snr,pertype)
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
Let \eqn{I(f_i)} and \eqn{S(f_i)} be respectively the periodogram of a given perturbed fractional Gaussian noise and the spectral density of perturbed fGn with Hurst exponent \code{eta} and the signal-to-noise ratio \code{SNR} where \eqn{f_i=2\pi *i/n} with \eqn{i=1,...,(n-1)/2}. Then, the value of log periodogram variance function is determined as

\eqn{\epsilon (f)\equiv \log \left ( \frac{I(f)}{S_X(f)} \right )+\gamma }.

See Percival and Walden (2000) for details. Some parts of this function were adopted from the S-PLUS codes originally developed by Jan Beran.
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

Percival and Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge University Press.

Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall.
}
\author{ Wonsang You }

\seealso{ \link{Qeta}, \link{Csum}, \link{QLCfun} }
\examples{
n<-1000; H<-0.7; SNR<-10
ts <- perturbFGN(n,H,type="WN",SNR=SNR)
yper<-per(ts)
ts.lpvar<-Lpvar(H,n,yper,snr=SNR,pertype="per")
}