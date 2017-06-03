\name{wFGN}
\alias{wFGN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Whittle estimator for a contaminated fGn }
\description{
 Whittle estimator of the Hurst parameter of a fractional Gaussian noise contaminated by additive outliers or noise.
}
\usage{
wFGN(x, istart=1, iend=length(x), nloop=1, init=c(0.55,0.01), 
		ndeps=c(1e-7,1e-2), noise=TRUE, pertype="per", minfun="qeta", 
		weights=c(1,1,0), cluster=FALSE, print.level=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a time series }
  \item{istart}{ the index of start. Default is \code{1} }
  \item{iend}{ the index of end. Default is the length of \code{x}. }
  \item{nloop}{ the number of sub-series. The Hurst parameter of each sub-series is estimated respectively. Default is \code{1}. }
  \item{init}{ the initial search values of Hurst parameter and SNR. Default is \code{c(0.55,0.01)}. }
  \item{ndeps}{ A vector of step sizes for optimization. Default is \code{c(1e-7,1e-2)}. }
  \item{noise}{ Enable the assumption of noise corruption. Default is \code{TRUE}. }
  \item{pertype}{ type of periodogram. Possible modes are \code{"per","taper"}. }
  \item{minfun}{ type of minimization function. Possible modes are \code{"qeta","lp","csum"}, and \code{"combi"}. }
  \item{weights}{ A vector of weights for each minimization function when \code{minfun="combi"}. }
  \item{cluster}{ A vector of machine names for parallel processing. For details, refer to the manual of package \code{genoud}. }
  \item{print.level}{ The logic value which controls the printing level. Four levels are available: 0 (minimal), 1 (normal), 2 (detailed), and 3 (debug). Default is 2. Refer to the manual of package \code{rgenoud} for details. }
}

\details{
The Hurst parameter of a pure or corrupted fractional Gaussian noise is estimated by a modified Whittle estimator. The most parts of this function were adopted from the S-PLUS codes originally developed by Jan Beran. See Beran (1994) for details.
}
\value{
 \item{thetavector}{ a vector of scaling coefficient and Hurst parameter }
 \item{Hlow}{ the lower limit of \eqn{95\%} confidence interval of Hurst parameter }
 \item{Hup}{ the upper limit of \eqn{95\%} confidence interval of Hurst parameter }
 \item{SNR}{ the signal-to-noise ratio }
 \item{fest}{ the spectral density of an estimated fGn }
}

\references{
Wonsang You (2010) Modified Whittle's Maximum Likelihood Estimator for Fractional Gaussian Noises Contaminated by Additive Noises, Technical Reports of the Leibniz Institute for Neurobiology, TR10015.

Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall.
}
\author{ Wonsang You }

\seealso{ \code{\link{circFGN}}, \code{\link{perturbFGN}}, \code{\link{wFGN.eval}} }

\examples{
n<-1000; H1<-0.2; H2<-0.7
ts1 <- circFGN(n,H1)
dat1 <- wFGN(ts1,noise=FALSE)
ts2 <- perturbFGN(n,H2,type="WN",SNR=10)
dat2 <- wFGN(ts2,pertype="per",minfun="combi")
}