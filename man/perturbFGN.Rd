\name{perturbFGN}
\alias{perturbFGN}

\title{ Simulation of a perturbed fGn }
\description{
Simulation of a sample path of a fractional Gaussian noise contaminated by outliers or an additive Gaussian noise 
}
\usage{
perturbFGN(n, H, C = 1, type = "no", p=.005, SNR=NULL, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ sample size }
  \item{H}{ Hurst parameter in (0,1) }
  \item{C}{ scaling coefficient. Default is \code{1}. }
  \item{type}{ type of perturbation. Possible modes are \code{"no","WN","AO"}.}
  \item{p}{ the ratio of perturbation by outliers. Default is 0.005. }
  \item{SNR}{ Signal to noise ratio }
  \item{plot}{ plot of perturbed fGn. Default is \code{FALSE}.}
}
\details{
The types of contamination are listed below. The most parts of this function were adopted from the package \code{dvfBm}.
\describe{
\item{"no"}{ no contamination}
\item{"AO"}{ contamination by additive outliers of 0.5\%}
\item{"WN"}{ contamination by i.i.d. standard Gaussian noise }
}
}
\value{
  a vector of length \code{n}.
}
\references{ 
S. Achard and J.-F. Coeurjolly (2009). Discrete variations of the fractional Brownian in the presence of outliers and an additive noise.
}
\author{ Wonsang You }

\seealso{ \code{\link{circFGN}}, \code{\link{wFGN}} }
\examples{
ts1<-perturbFGN(1000,0.2,type="no",plot=TRUE)
ts2<-perturbFGN(1000,0.7,type="AO",SNR=20,plot=TRUE)
ts3<-perturbFGN(10000,0.3,type="WN",SNR=10,plot=TRUE)
}