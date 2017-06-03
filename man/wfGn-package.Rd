\name{wfGn-package}
\alias{wfGn-package}
\docType{package}
\title{
Estimation of the Hurst Parameter for Contaminated Fractional Gaussian Noises
}
\description{Estimates the Hurst parameter of a contaminated fractional Gaussian noises by using the modified Whittle estimator.
}
\details{
\tabular{ll}{
Package: \tab wfGn\cr
Type: \tab Package\cr
Version: \tab 1.2-1\cr
Date: \tab 2010-10-29\cr
License: \tab GPL (>=2.0)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Wonsang You and Jan Beran
}
\references{
Wonsang You (2010) Modified Whittle's Maximum Likelihood Estimator for Fractional Gaussian Noises Contaminated by Additive Noises, Technical Reports of the Leibniz Institute for Neurobiology, TR10015.

Percival and Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge University Press.

Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall. 

S. Achard and J.-F. Coeurjolly (2009). Discrete variations of the fractional Brownian in the presence of outliers and an additive noise.
}

\keyword{ package }

\examples{
n<-10000;H<-0.2
z<-perturbFGN(n,H,type="WN",SNR=-10,plot=TRUE)
wFGN(z)
}