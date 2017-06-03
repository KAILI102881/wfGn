\name{QLCfun}
\alias{QLCfun}

\title{ Combinational Minimization Function }
\description{
Combinational function of Qeta, Cumulative Qeta Difference, and Log Periodogram Variance Function to be minimized to estimate the parameters of a fractional Gaussian noise.
}
\usage{
QLCfun(eta,n,yper,snr,pertype,weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta}{ a positive value of the Hurst exponent which is less than 1. }
  \item{n}{ the number of time points. }
  \item{yper}{ a vector of periodogram with length of the largest integer less than \code{(n-1)/2} }
  \item{snr}{ the signal-to-noise ratio. }
  \item{pertype}{ the type of periodogram. Possible modes are \code{"per","taper"}. }
  \item{weights}{ a vector of length 3 which indicates the weights of Qeta, Cumulative Qeta Difference, and Log Periodogram Variance Functions }
}
\details{
Let \eqn{\widetilde{L}_1(y;\theta)}, \eqn{\widetilde{L}_2(y;\theta)}, \eqn{\widetilde{L}_3(y;\theta)} be Qeta, Cumulative Qeta Difference, and Log Periodogram Variance Functions respectively. Also, let \eqn{p=\left ( p_1,p_2,p_3 \right )} be a weight vector, and let \eqn{\widetilde{L}\left ( y;\theta \right )} be a combinational minimization function. Then, the value of the combinational function is determined as

\eqn{\widetilde{L}\left ( y;\theta \right )=p_1 \times \widetilde{L}_1\left ( y;\theta \right )+p_2 \times \widetilde{L}_2\left ( y;\theta \right )+p_3 \times \widetilde{L}_3\left ( y;\theta \right )}.

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

Percival and Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge University Press.

Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall.
}
\author{ Wonsang You }

\seealso{ \code{\link{Qeta}}, \code{\link{Csum}}, \code{\link{Lpvar}} }
\examples{
n<-1000; H<-0.7; SNR<-10
ts <- perturbFGN(n,H,type="WN",SNR=SNR)
yper<-per(ts)
ts.qlc<-QLCfun(H,n,yper,snr=SNR,pertype="per")
}