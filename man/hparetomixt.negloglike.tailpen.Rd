\name{hparetomixt.negloglike.tailpen}
\alias{hparetomixt.negloglike.tailpen}
\alias{hparetomixt.fit.tailpen}
\alias{hparetomixt.cvtrain.tailpen}
\title{
Maximum Likelihood Estimation for a Mixture of Hybrid Paretos with Tail Penalty
}
\description{
In order to drive the tail index estimation, a penalty is introduced in
the log-likelihood. The goal of the penalty is to include a priori
information which in our case is that only a few mixture components have
a heavy tail index which should approximate the tail of the underlying
distribution while most other mixture components have a light tail and
aim at modelling the central part of the underlying distribution.
}
\usage{
hparetomixt.negloglike.tailpen(params, lambda, w, beta, mu, sigma, x)
hparetomixt.fit.tailpen(params, lambda, w, beta, mu, sigma, x, ...)
hparetomixt.cvtrain.tailpen(m, lambda, w, beta, mu, sigma, x, nfold=5, nstart=1, ...)
}
\arguments{
  \item{params}{matrix of dimension 4 by m, where m is the number of
    components, each column of the matrix contains the mixture
    parameters of one component (pi, xi, mu, sigma)
}
  \item{m}{number of mixture components}

  \item{lambda}{penalty parameter which controls the trade-off between
    the penalty and the negative log-likelihood, takes on positive
    values
    }
  \item{w}{penalty parameter in [0,1] which is the proportion of
    components with light tails, 1-\code{w} being the proportion of
    components with heavy tails
}
  \item{beta}{positive penalty parameter which indicates the importance
    of the light tail components (it is the parameter of an exponential
    which represents the prior over the light tail components)

  }

  \item{mu}{penalty parameter in (0,1) which represents the a priori
    value for the heavy tail index of the underlying distribution
}
  \item{sigma}{positive penalty parameter which controls the spread
  around the a priori value for the heavy tail index of the underlying
  distribution
 }
  \item{x}{a vector of length n of observations assumed to be sampled from a
    mixture of hybrid Paretos
  }
  \item{nfold}{number of fold for cross-validation estimate, default is 5}
  \item{nstart}{number of re-starts for the optimizer \code{nlm} with
    different initial parameters, default is 1}
  \item{\dots}{
optional arguments for \code{nlm}
}
}
\details{
  The penalty term is given by the logarithm of the following
    two-component mixture, as a function of a tail index parameter xi :
  w beta exp(-beta xi) + (1-w) exp(-(xi-mu)^2/(2 sigma^2))/(sqrt(2 pi)
  sigma)
  where the first term is the prior on the light tail component and the
    second term is the prior on the heavy tail component.
}
\value{
\code{hparetomixt.negloglike.tailpen} returns a single value (the negative log-likelihood for
given parameters and sample) and a vector, the gradient, which is passed as an attribute,
while \code{hparetomixt.fit.tailpen} returns a 4 by m matrix of MLE for the
hybrid Pareto mixture parameters and \code{hparetomixt.cvtrain.tailpen}
    returns a cross-validation estimate of the out-of-sample negative
    log-likelihood for the given model (number of components and penalty parameters)
}
\references{
 Carreau, J.,Naveau, P. and Sauquet, E. (2009), A statistical
    rainfall-runoff mixture model with heavy-tailed components, 45,
    Water Resources Research
}
\author{
Julie Carreau
}

\seealso{
\code{\link{hparetomixt.init}}, \code{\link{hparetomixt.negloglike}}
}
\examples{
r <- rfrechet(500,loc=5,scale=5,shape=5)
m <- 2
param.init <- hparetomixt.init(m,r)
hparetomixt.negloglike.tailpen(param.init,10,0.5,20,0.1,0.2,r)
hparetomixt.fit.tailpen(param.init,10,0.5,20,0.1,0.2,r)
hparetomixt.cvtrain.tailpen(2,10,0.5,20,0.1,0.2,r)
}
