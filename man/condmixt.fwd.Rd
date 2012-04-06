\name{condmixt.fwd}
\alias{condmixt.fwd}
\alias{condhparetomixt.fwd}
\alias{condhparetomixt.dirac.fwd}
\alias{condgaussmixt.fwd}
\alias{condgaussmixt.dirac.fwd}
\alias{condbergamixt.fwd}
\title{
Forward pass in neural network conditional mixtures
}
\description{
A forward pass means that given explanatory variables \code{x}, the
neural network computes the corresponding values of the mixture
parameters.
}
\usage{
condhparetomixt.fwd(theta, h, m, x)
condhparetomixt.dirac.fwd(theta, h, m, x)
condgaussmixt.fwd(theta, h, m, x)
condgaussmixt.dirac.fwd(theta, h, m, x)
condbergamixt.fwd(theta,h,x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{Vector of neural network parameters
}
  \item{h}{Number of hidden units
}
  \item{m}{Number of components
}
\item{x}{
  Matrix of explanatory (independent) variables of dimension d x n, d is
the number of variables and n is the number of examples (patterns)
}
}
\details{
  \code{condhparetomixt} indicates a mixture with hybrid Pareto
components,
\code{condgaussmixt} for Gaussian components,
\code{condbergam} for a Bernoulli-Gamma two component mixture,
\code{dirac} indicates that a discrete dirac component is included in
the mixture
The forward pass for Log-Normal conditional mixture is the same one as
for Gaussian conditional mixture. Therefore, \code{condgaussmixt.fwd} can be
used for the forward pass of a conditional mixture with Log-Normal components.
}
\value{
  A matrix of mixture parameters corresponding to the values in \code{x} :
  - for \code{condhparetomixt.fwd}, each component requires four
  parameters (pi, xi, mu, sigma) and the parameter matrix has dimensions
  \code{m} x 4 x n
  - for \code{condhparetomixt.dirac.fwd}, there is an additional
  parameter for the probability of the dirac at zero so that the mixture
  parameters are stored in a (4\code{m}+1) x n matrix
  - for \code{condgaussmixt.fwd}, each component requires three
  parameters (pi, mu, sigma) and the parameter matrix has dimensions
  \code{m} x 3 x n
  - for \code{condgaussmixt.dirac.fwd}, there is an additional
  parameter for the probability of the dirac at zero so that the mixture
  parameters are stored in a (3\code{m}+1) x n matrix
  - for \code{condbergamixt.fwd}, there are three parameters, the
  probability of the dirac at zero and two parameters for the Gamma distribution
  }
  \references{
    Bishop, C. (1995), Neural Networks for Pattern Recognition, Oxford

    Carreau J. and Vrac, M. (2011) Stochastic Downscaling of Precipitation
    with Neural Network Conditional Mixture Models, 47, Water Resources Research

    Williams, M.P. (1998) Modelling Seasonality and Trends in Daily Rainfall
    Data, 10, Advances in Neural Information and Processing Systems
}
\author{
Julie Carreau
}

\seealso{
\code{\link{condmixt}} \code{\link{condmixt.init}}, \code{\link{condmixt.nll}}
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

params.mixt <- condhparetomixt.fwd(thetainit,h,m,t(x)) # compute mixture parameters

}
