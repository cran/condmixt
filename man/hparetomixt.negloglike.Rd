\name{hparetomixt.negloglike}
\alias{hparetomixt.negloglike}
\alias{hparetomixt.fit}
\alias{hparetomixt.cvtrain}
\title{
Maximum Likelihood Estimation for a Mixture of Hybrid Paretos
}
\description{
Negative log-likelihood and gradient (\code{hparetomixt.negloglike}),
MLE of a hybrid Pareto distribution parameters
(\code{hparetomixt.fit}) and out-of-sample negative log-likelihood
estimation for a given number of components with nfold cross-validation
(\code{hparetomixt.cvtrain}).

\code{hparetomixt.fit} applies the optimizer \code{nlm} to minimize the
negative log-likelihood based on some starting values for the hybrid
Pareto parameters.
}
\usage{
hparetomixt.negloglike(params, x)
hparetomixt.fit(params, x, ...)
hparetomixt.cvtrain(m, x, nfold=5, nstart=1, ...)
}
\arguments{
  \item{params}{matrix of dimension 4 by m, where m is the number of
    components, each column of the matrix contains the mixture
    parameters of one component (pi, xi, mu, sigma)
}
  \item{x}{
a vector of length n of observations assumed to be sampled from a
mixture of hybrid Paretos
}
\item{m}{number of mixture components}
\item{nfold}{number of fold for cross-validation estimate, default is 5}
\item{nstart}{number of re-starts for the optimizer \code{nlm} with
  different initial parameters, default is 1}
\item{\dots}{
optional arguments for \code{nlm}
}
}
\value{
\code{hparetomixt.negloglike} returns a single value (the negative log-likelihood for
given parameters and sample) and a vector, the gradient, which is passed as an attribute,
while \code{hparetomixt.fit} returns a 4 by m matrix of MLE for the
hybrid Pareto mixture parameters and \code{hparetomixt.cvtrain}
    returns a cross-validation estimate of the out-of-sample negative
    log-likelihood for a selected number of components
}
\references{
Carreau, J. and Bengio, Y. (2009), A Hybrid Pareto Model for
Asymmetric Fat-tailed Data: the Univariate Case, 12, Extremes
}
\author{
Julie Carreau
}
\seealso{
\code{\link{hparetomixt.init}}, \code{\link{hparetomixt.negloglike.tailpen}}
}
\examples{
r <- rfrechet(500,loc=5,scale=5,shape=5)
m <- 2
param.init <- hparetomixt.init(m,r)
hparetomixt.negloglike(param.init,r)
hparetomixt.fit(param.init,r)
}
