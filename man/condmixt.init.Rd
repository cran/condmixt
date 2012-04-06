\name{condmixt.init}
\alias{condmixt.init}
\alias{condhparetomixt.init}
\alias{condhparetomixt.dirac.init}
\alias{condgaussmixt.init}
\alias{condgaussmixt.dirac.init}
\alias{condlognormixt.init}
\alias{condlognormixt.dirac.init}
\alias{condbergamixt.init}
\title{
Conditional mixture parameter initial values
}
\description{
Neural network weights are randomly initialized uniformly over the range
[-0.9/sqrt(k),0.9/sqrt(k)] where k is the number of inputs to the
neuron. This ensures that the hidden units will not be saturated and
that training should proceed properly. In addition, if the dependent data
Y is provided, the biases will be initialized according to the initial
parameters of an unconditional mixture computed on the dependent data.
}
\usage{
condhparetomixt.init(d, h, m, y = NULL)
condhparetomixt.dirac.init(d, h, m, y = NULL)
condgaussmixt.init(d,h,m,y=NULL)
condgaussmixt.dirac.init(d,h,m,y=NULL)
condlognormixt.init(d,h,m,y=NULL)
condlognormixt.dirac.init(d,h,m,y=NULL)
condbergamixt.init(d,h,y=NULL)
}
\arguments{
  \item{d}{
dimension of input x to neural network
}
  \item{h}{
number of hidden unit
}
  \item{m}{
number of components
}
  \item{y}{
optional, dependent one-dimensional data
}
}
\details{
If the argument \code{y} is provided, an unconditional mixture with the
same type of components will be initialized on \code{y}. These initial
unconditional parameters are then used to give more appropriate initial
values to the biases of the neural network.
}
\value{
A vector of neural network parameters for the given number of hidden
units and number of components and specific conditional mixture
formulation : hybrid Pareto components (condhparetomixt.init), hybrid
Pareto components + discrete dirac component at zero
(condhparetomixt.dirac.init), Gaussian components (condgaussmixt.init),
Gaussian components + discrete dirac component at zero
(condgaussmixt.dirac.init), Log-Normal components (condlognormixt.init),
Log-Normal components + discrete dirac component at zero
(condlognormixt.dirac.init) and the Bernoulli-Gamma two-component
mixture (condbergamixt.init).
}
\references{
Carreau, J. and Bengio, Y. (2009), A Hybrid {P}areto Mixture for
Conditional Asymmetric Fat-Tailed Distributions, 20, IEEE Transactions
on Neural Networks

Nabney, I. (2002) NetLab : Algorithms for Pattern Recognition, Springer

Williams, M.P. (1998) Modelling Seasonality and Trends in Daily Rainfall
Data, 10, Advances in Neural Information and Processing Systems
}
\author{
Julie Carreau
}
\seealso{
\code{\link{hparetomixt.init}}, \code{\link{gaussmixt.init}}
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
}

