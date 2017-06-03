\name{wFGN.eval}
\alias{wFGN.eval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Evaluation of the modified Whittle estimator for a contaminated fGn }
\description{
 Statistical performance evaluation of the modified Whittle estimator of the Hurst parameter of a fractional Gaussian noise contaminated by additive outliers or noise.
}
\usage{
wFGN.eval(H=NULL, n=1000, m=100, type="no", SNR=NULL, 
			ndeps=c(1e-7,1e-2), noise=TRUE, pertype="per", 
			minfun="qeta", weights=c(1,1,0), cluster=FALSE, 
			plot=TRUE, sav=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H}{ Hurst parameter. If \code{H=NULL}, \eqn{H=0.1,0.2,...,0.9} are tested. Default is \code{NULL}. }
  \item{n}{ sample size. Default is \code{1000}. }
  \item{m}{ the number of repetitions. Default is \code{100}. }
  \item{type}{ type of perturbation. Possible modes are \code{"no","WN","AO"}. Refer to \link{wFGN} for details. }
  \item{SNR}{ Signal to noise ratio. }
  \item{ndeps}{ A vector of step sizes for optimization. Default is \code{c(1e-7,1e-2)}. }
  \item{noise}{ Enable the assumption of noise corruption. Default is \code{TRUE}. }
  \item{pertype}{ type of periodogram. Possible modes are \code{"per","taper"}. }
  \item{minfun}{ type of minimization function. Possible modes are \code{"qeta","lp","csum","combi"}. }
  \item{weights}{ A vector of weights for each minimization function when \code{minfun="combi"}. }
  \item{cluster}{ A vector of machine names for parallel processing. For details, refer to the manual of package \code{genoud}. }
  \item{plot}{ a boxplot of parameter estimation. Default is \code{TRUE}. }
  \item{sav}{ Enable sample plots. Default is \code{FALSE}. }
}

\details{
The Hurst parameter of a fractional Gaussian noise is estimated by the modified Whittle estimator. This function evaluates the consistency of the Whittle estimator by several repetitions.
}
\value{
\item{Hdata}{ a \eqn{m \times 1} or \eqn{m \times 9} matrix of Hurst parameter estimates for fGn with different Hurst parameters}
\item{Hstat}{ a  \eqn{4 \times 1} or \eqn{4 \times 9} matrix with a sample Hurst parameter, mean, standard deviation, and mean squared error (MSE) of Hurst parameter estimates }
\item{SNRdata}{ a \eqn{m \times 1} or \eqn{m \times 9} matrix of SNR estimates for fGn with different Hurst parameters}
\item{SNRstat}{ a  \eqn{4 \times 1} or \eqn{4 \times 9} matrix with a sample SNR, mean, standard deviation, and mean squared error (MSE) of SNR estimates }
\item{Theta}{ a \eqn{m \times 1} or \eqn{m \times 9} matrix of scaling coefficient estimates }
}

\references{
Wonsang You (2010) Modified Whittle's Maximum Likelihood Estimator for Fractional Gaussian Noises Contaminated by Additive Noises, Technical Reports of the Leibniz Institute for Neurobiology, TR10015.
}
\author{ Wonsang You }

\seealso{ \code{\link{circFGN}}, \code{\link{perturbFGN}}, \code{\link{wFGN}} }

\examples{
dat <- wFGN.eval(H=0.2, n=10000, m=10, type="WN", SNR=-20, plot=TRUE)
}