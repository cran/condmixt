\name{hparetomixt.init}
\alias{hparetomixt.init}
\title{
Provides initial values for the parameters of a mixture of hybrid
Paretos based on a sample.
}
\description{
Initial values for the parameters of a mixture of hybrid Paretos are
provided by applying the following steps :
1) clustering the sample into as many clusters as there are mixture
components
2) estimating the hybrid Pareto parameters for each component on the data
from each cluster with the moment-like estimators, see
\code{\link{hpareto.mme}}
}
\usage{
hparetomixt.init(m, x, iter.max = 20, nstart = 10)
}
\arguments{
  \item{m}{number of mixture components}
  \item{x}{data sample from which the initial parameters are computed}
  \item{iter.max}{maximum number of iteration for kmeans clustering,
    default is 20, see \code{\link{kmeans}}}
  \item{nstart}{number of random cluster centers chosen (default is 10), see
    \code{\link{kmeans}}}
}
\value{
a matrix of dimension 4 x \code{m} which stores the 4 parameters (pi,
xi, mu, sigma) of each of the \code{m} components.
}
\references{
Carreau, J. and Bengio, Y. (2009), A Hybrid Pareto Model for
Asymmetric Fat-tailed Data: the Univariate Case, 12, Extremes
}
\author{
Julie Carreau
}
\seealso{
  \code{\link{kmeans}}, \code{\link{hpareto.mme}}
}
\examples{
r <- rfrechet(500,loc=5,scale=5,shape=5)
m <- 2
param.init <- hparetomixt.init(m,r)
}
