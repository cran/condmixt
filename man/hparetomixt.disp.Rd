\name{hparetomixt.disp}
\alias{hparetomixt.disp}
\title{
Display the Hybrid Pareto Mixture Parameters}
\description{
 Display the Hybrid Pareto Mixture Parameters in Latex Style Tabular form
}
\usage{
hparetomixt.disp(params)
}
\arguments{
  \item{params}{matrix of mixture parameters of dimension 4 x m, where m
    is the number of components, so that each column contains the
    mixture parameters (pi, xi, mu, sigma) related to one component
}
}
\value{
Display on the R consol the hybrid Pareto mixture parameters with syntax
suitable to make a Latex tabular out of it.
}
\author{
Julie Carreau
}
\seealso{
\code{\link{hparetomixt}}, \code{\link{hparetomixt.negloglike}}
}
\examples{
r <- rfrechet(500,loc=5,scale=5,shape=5)
m <- 2
param.init <- hparetomixt.init(m,r)
hparetomixt.disp(param.init)
}
