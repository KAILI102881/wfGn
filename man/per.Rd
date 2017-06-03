\name{per}
\alias{per}

\title{ Periodogram Estimator }
\description{
Estimation of periodogram of a given time series.
}
\usage{
per(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ a vector which indicates a time series }
}
%\details{
%}
\value{
  a vector of periodogram with length of \eqn{1+ floor(n/2)} where \eqn{n} is the length of \eqn{z}. This function is originated from Jan Beran. See Beran (1994) for details.
}
\references{ 
Jan Beran (1994) Statistics for Long-Memory Processes, Chapman & Hall.
}
\author{ Jan Beran }

%\seealso{}
\examples{
n<-1000; H<-0.7;
ts <- circFGN(n,H)
yper<-per(ts)
}