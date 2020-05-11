\name{condmixt.fit}
\alias{condmixt.fit}
\alias{condhparetomixt.fit}
\alias{condhparetomixt.fit.tailpen}
\alias{condhparetomixt.dirac.fit}
\alias{condhparetomixt.dirac.fit.tailpen}
\alias{condgaussmixt.fit}
\alias{condgaussmixt.dirac.fit}
\alias{condlognormixt.fit}
\alias{condlognormixt.dirac.fit}
\alias{condbergamixt.fit}
\title{
Maximum likelihood estimation for conditional mixture parameters
}
\description{
Performs maximum likelihood estimation of conditional mixture parameters
given starting values and the numbers of hidden units and of components, the
type of components and the presence of a discrete dirac component on a
given data set.
}
\usage{
condhparetomixt.fit(theta, h, m, x, y, ...)
condhparetomixt.fit.tailpen(theta,h,m,x,y,lambda=0,w=0.2,beta=50,
                            mu=0.2,sigma=0.2, ...)
condhparetomixt.dirac.fit.tailpen(theta,h,m,x,y,lambda=0,w=0.2,beta=50,
                                  mu=0.2,sigma=0.2, ...)
condhparetomixt.dirac.fit(theta,h,m,x,y, ...)
condgaussmixt.fit(theta,h,m,x,y, ...)
condgaussmixt.dirac.fit(theta,h,m,x,y, ...)
condlognormixt.fit(theta,h,m,x,y, ...)
condlognormixt.dirac.fit(theta,h,m,x,y, ...)
condbergamixt.fit(theta,h,x,y, ...)
}
\arguments{
  \item{theta}{
Vector of neural network parameters
}
  \item{h}{
Number of hidden units
}
  \item{m}{
Number of components
}
  \item{x}{
Matrix of explanatory (independent) variables of dimension d x n, d is
the number of variables and n is the number of examples (patterns)
}
  \item{y}{
Vector of n dependent variables
}
 \item{lambda}{penalty parameter which controls the trade-off between
    the penalty and the negative log-likelihood, takes on positive
    values. If zero, no penalty
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
\item{...}{optional arguments for the optimizer \code{nlm}}
}
\details{
\code{condhparetomixt} indicates a mixture with hybrid Pareto
components,
\code{condgaussmixt} for Gaussian components,
\code{condlognormixt} for Log-Normal components,
\code{condbergam} for a Bernoulli-Gamma two component mixture,
\code{tailpen} indicates that a penalty is added to the log-likelihood
to guide the tail index parameter estimation,
\code{dirac} indicates that a discrete dirac component at zero is included in
the mixture

In order to drive the tail index estimation, a penalty is introduced in
the log-likelihood. The goal of the penalty is to include a priori
information which in our case is that only a few mixture components have
a heavy tail index which should approximate the tail of the underlying
distribution while most other mixture components have a light tail and
aim at modelling the central part of the underlying distribution.

The penalty term is given by the logarithm of the following
two-component mixture, as a function of a tail index parameter xi :

w beta exp(-beta xi) + (1-w) exp(-(xi-mu)^2/(2 sigma^2))/(sqrt(2 pi) sigma)

where the first term is the prior on the light tail component and the
second term is the prior on the heavy tail component.
}
\value{Returns a vector of neural network parameters corresponding to
  the maximum likelihood estimation.
}
\references{
Bishop, C. (1995), Neural Networks for Pattern Recognition, Oxford

Carreau, J. and Bengio, Y. (2009), A Hybrid {P}areto Mixture for
Conditional Asymmetric Fat-Tailed Distributions, 20, IEEE Transactions
on Neural Networks
}
\author{
Julie Carreau
}
\seealso{
\code{\link{condmixt.nll}}, \code{\link{condmixt.init}}
}
\examples{
n <- 200
x <- runif(n) # x is a random uniform variate
# y depends on x through the parameters of the Frechet distribution
y <- rfrechet(n,loc = 3*x+1,scale = 0.5*x+0.001,shape=x+1)

plot(x,y,pch=22)
# 0.99 quantile of the generative distribution
qgen <- qfrechet(0.99,loc = 3*x+1,scale = 0.5*x+0.001,shape=x+1)
points(x,qgen,pch="*",col="orange")

h <- 2 # number of hidden units
m <- 4 # number of components

# initialize a conditional mixture with hybrid Pareto components
thetainit <- condhparetomixt.init(1,h,m,y)

# MLE for initial neural network parameters
condhparetomixt.fit(thetainit,h,m,t(x),y,iterlim=100)
}
