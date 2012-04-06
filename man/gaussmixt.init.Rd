\name{gaussmixt.init}
\alias{gaussmixt.init}
\title{
  Provides initial values for the parameters of a mixture of Gaussians
  based on a sample.
}
\description{
Initial values for the parameters of a mixture of Gaussians are
provided by applying the following steps :
1) clustering the sample into as many clusters as there are mixture
components
2) the initial means and standard deviations of each component are taken
as the cluster centers and median absolute deviation respectively
computed on each component
}
\usage{
gaussmixt.init(m, x, iter.max = 20, nstart = 10)
}
\arguments{
  \item{m}{number of mixture components
}
  \item{x}{data sample from which the initial parameters are computed
  }
  \item{iter.max}{maximum number of iteration for kmeans clustering,
    default is 20, see \code{\link{kmeans}}}
  \item{nstart}{number of random cluster centers chosen (default is 10), see
    \code{\link{kmeans}}}
}

\value{
  a matrix of dimension 3 x \code{m} which stores the 3 parameters (pi,
  mu, sigma) of each of the \code{m} components.
}
\references{
McLachlan, G. and Peel, D. (2000), Finite Mixture Models, Wiley series
  in probability and statistics
}
\author{
Julie Carreau
}
\seealso{
\code{\link{kmeans}},\code{\link{hparetomixt.init}}
}
\examples{
r <- rfrechet(500,loc=5,scale=5,shape=10)
m <- 2
param.init <- gaussmixt.init(2,r)
}
