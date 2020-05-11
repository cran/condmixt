\name{hpareto.negloglike}
\alias{hpareto.negloglike}
\alias{hpareto.fit}
\title{
Hybrid Pareto Maximum Likelihood Estimation
}
\description{
Negative log-likelihood and gradient (\code{hpareto.negloglike}) and MLE of a hybrid Pareto
distribution parameters (\code{hpareto.fit}). \code{hpareto.fit} applies
the optimizer \code{nlm} to minimize the negative log-likelihood based on
some starting values for the hybrid Pareto parameters.
}
\usage{
hpareto.negloglike(params, x)
hpareto.fit(params, x, ...)
}
\arguments{
  \item{params}{
vector of length 3 of hybrid Pareto parameters c(xi, mu, sigma)
}
  \item{x}{
a vector of length n of observations assumed to be sampled from a hybrid Pareto distribution
}
  \item{\dots}{
optional arguments for \code{nlm} called by \code{hpareto.fit}
}
}

\value{
\code{hpareto.negloglike} returns a single value (the negative log-likelihood for
given parameters and sample) and a vector, the gradient, which is passed as an attribute,
while \code{hpareto.fit} returns a vector of three parameters, the MLE
of the parameters of the hybrid Pareto distribution
}
\references{
Carreau, J. and Bengio, Y. (2009), A Hybrid Pareto Model for
Asymmetric Fat-tailed Data: the Univariate Case, 12, Extremes
}
\author{
Julie Carreau
}

\seealso{
\code{\link{hpareto.mme}}, \code{\link{dhpareto}},\code{\link{phpareto}} and \code{\link{rhpareto}}
}
\examples{
r <- rhpareto(500,0.2,trunc=FALSE)
params.init <- hpareto.mme(r)
hpareto.negloglike(params.init,r)
hpareto.fit(params.init,r)
}