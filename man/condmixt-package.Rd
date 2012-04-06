\name{condmixt-package}
\alias{condmixt-package}
\alias{condmixt}
\docType{package}
\title{
Conditional Density Estimation with Neural Network Conditional Mixtures
}
\description{
Neural network conditional mixtures are mixtures whose parameters
depends on explanatory variables through a neural network. In other
words, for a given set of values of explanatory variables, a neural
network will compute the corresponding mixture parameters. Thus, the
parameters of the mixture and hence the shape of its density is modified
by the values of the explanatory variables. There are two
hyper-parameters, the number of hidden units of the neural network and
the number of components of the mixture, which control the complexity of
the conditional mixtures. In this package, a number of types of
conditional mixtures are implemented, which differ in their type of
components and which allow the possibility to include a discrete
component in the mixture. 
}
\details{
\tabular{ll}{
Package: \tab condmixt\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-04-06\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Functions for conditional mixtures with
- hybrid Pareto components start with \code{condhparetomixt} 
- Gaussian components start with \code{condgaussmixt}
- Log-Normal components start with \code{condlognormixt}

If a discrete dirac component at zero is added, for instance to model the
presence of zeros for "no rain" events in a rainfall time-series, then
functions related to such mixtures start with something like
\code{condhparetomixt.dirac} depending on the type of components for the
continuous part.

One special mixture, which has only two-components (a discrete, dirac or
Bernoulli, component and a Gamma component) is implemented under the
name \code{condbergamixt}. It has been introduced in Williams(1998).

Finally, hybrid Pareto components have a tail index parameter which
might be difficult to estimate, specially in the conditional,
non-stationnary case. To alleviate this issue, a penalty term might be
added to the log-likelihood in order to guide maximum likelihood
estimation of the tail indexes. The functions which employ a tail
penalty term have a name ending with \code{tailpen}. 

The goal of the tail penalty for hybrid Pareto mixtures is to include a priori
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
\author{
Julie Carreau

Maintainer: Julie Carreau <julie.carreau@univ-montp2.fr>
}
\references{
  Bishop, C. (1995), Neural Networks for Pattern Recognition, Oxford

  Carreau J. and Vrac, M. (2011) Stochastic Downscaling of Precipitation
  with Neural Network Conditional Mixture Models, 47, Water Resources Research

Carreau, J. and Bengio, Y. (2009), A Hybrid Pareto Model for Asymmetric
Fat-tailed Data: the Univariate Case, 12, Extremes

Carreau, J. and Bengio, Y. (2009), A Hybrid {P}areto Mixture for
Conditional Asymmetric Fat-Tailed Distributions, 20, IEEE Transactions
on Neural Networks

Williams, M.P. (1998) Modelling Seasonality and Trends in Daily Rainfall
Data, 10, Advances in Neural Information and Processing Systems
}
\seealso{
\code{\link{condmixt.init}}, \code{\link{condmixt.nll}}, \code{\link{condmixt.fit}}
}

