\name{hillest}
\alias{hillest}
\title{
Hill Estimator
}
\description{
Hill estimator of the tail index. This estimator assumes the tail index
to be positive. The threshold used is the k+1 order statistic.
}
\usage{
hillest(data, k)
}
\arguments{
  \item{data}{ vector of sample from a distribution for which the tail
    index is to be estimated
}
  \item{k}{define the number of order statistics used for the estimation
}
}
\value{Hill estimator of the tail index parameter.
}
\references{
Embrechts, P., Kluppelberg, C. and Mikosch, T. (1997), Modelling
Extremal Events, Applications of Mathematics, Stochastic Modelling and
Applied Probability, Springer
}
\author{
Julie Carreau
}
\seealso{
\code{\link{gpd.mme}}, \code{\link{hpareto.mme}}
}
\examples{
r<-rhpareto(5000,0.2,-50,1,trunc=TRUE)  # tail index is equal to 0.2
 # Hill estimator of the tail index with the 100 largest observations
hillest(r,10)
}
