\name{gpd.mme}
\alias{gpd.mme}
\alias{hpareto.mme}
\title{
Moment Estimator for the Generalized and the Hybrid Pareto Distribution
}
\description{
Moment estimators for the generalized Pareto distribution and parameter estimators
based on two quantiles plus a tail index estimator for the hybrid Pareto
distribution.
}
\usage{
gpd.mme(x)
hpareto.mme(x, xi0=c(),p=0.99)
}
\arguments{
  \item{x}{ vector of sample for which the parameters will be estimated}
  \item{xi0}{ an optional a priori value for the tail index parameter}
  \item{p}{ the percentage of largest observations used for tail index estimation  }
}
\details{
For \code{hpareto.mme}, the tail index is assumed to be positive. In
case one has some prior information on the value of the tail index
parameter, it is possible to provide this value as an argument to the
function \code{hpareto.mme}. The two other parameters mu and sigma will
be estimated based on that tail index estimate and two quantiles of the
hybrid Pareto distribution. If the tail index parameter is not provided
as input, it will be estimated with the Hill estimator using data above
the \code{p}-quantile. By default, \code{p}=0.99 but this might be inappropriate
depending on the sample. Since the tail index estimate is very
sensitive, it is recommended to tune carefully the \code{p} argument.
}
\value{
A vector of parameter estimates.
}
\references{
Hosking, J. R. M. and Wallis, J. R. (1987), Parameter and quantile
estimation for the Generalized Pareto distribution, 29, Technometrics
Carreau, J. and Bengio, Y. (2009), A Hybrid Pareto Model for
Asymmetric Fat-tailed Data: the Univariate Case, 12, Extremes
}
\author{
Julie Carreau
}
\examples{
r<-rhpareto(1000,0.1,0,1,trunc=FALSE)
hpareto.mme(r,p=0.991)
}
